"""
analysis_plot
=============

One cairo plot widget for the Trajectory Analysis window (matplotlib is not
an EasyHybrid dependency). Four kinds of content:

    plot.set_lines(series, xlabel, ylabel, title)       curves, points or bars
    plot.set_heatmap(matrix, ..., diverging=True)        matrix with a colour bar
    plot.set_categories(codes, colors, names, ...)       category map (DSSP)
    plot.set_scatter(x, y, c, ...)                       points coloured by c

series is a list of dicts: {"x": [...], "y": [...], "label": str,
"color": (r, g, b), "style": "line" | "points" | "bars", "dash": bool,
"colors": [(r, g, b), ...] (bars only, one per bar)}.

The value under the mouse is shown in a corner; on_click(x, y, info) is
called on a left click (the window uses it to jump to a frame or to select
residues in the 3D view). export(path) writes PNG or SVG.
"""

import math

import cairo
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk

PALETTE = [(0.12, 0.47, 0.71), (0.85, 0.37, 0.01), (0.17, 0.63, 0.17), (0.84, 0.15, 0.16),
           (0.58, 0.40, 0.74), (0.55, 0.34, 0.29), (0.89, 0.47, 0.76), (0.50, 0.50, 0.50),
           (0.74, 0.74, 0.13), (0.09, 0.75, 0.81)]


def nice_step(span, target=6):
    if span <= 0:
        return 1.0
    raw = span / target
    exp = math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if m * 10 ** exp >= raw:
            return m * 10 ** exp
    return 10 ** (exp + 1)


def _ticks(lo, hi, target=6):
    step = nice_step(hi - lo, target)
    if step <= 1e-12 * max(1.0, abs(lo), abs(hi)):
        return [lo, hi]
    start = math.ceil(lo / step - 1e-9) * step
    out = []
    v = start
    while v <= hi + 1e-9 * max(1.0, abs(hi)) and len(out) < 1000:
        out.append(round(v, 10))
        v += step
    return out


def _fmt(v):
    if v == 0:
        return "0"
    a = abs(v)
    if a >= 1e5 or a < 1e-3:
        return "{:.2e}".format(v)
    return "{:g}".format(round(v, 6))


def colormap(t, diverging=False):
    """t in [0, 1]. Sequential: white -> blue -> dark; diverging: blue - white - red."""
    t = min(1.0, max(0.0, t))
    if diverging:
        if t < 0.5:
            s = t / 0.5
            return (0.13 + 0.87 * s, 0.30 + 0.70 * s, 0.68 + 0.32 * s)
        s = (t - 0.5) / 0.5
        return (1.0 - 0.30 * s, 1.0 - 0.85 * s, 1.0 - 0.85 * s)
    stops = [(0.0, (1.0, 1.0, 0.85)), (0.35, (0.50, 0.80, 0.73)), (0.7, (0.17, 0.50, 0.72)), (1.0, (0.03, 0.15, 0.40))]
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t <= t1:
            s = (t - t0) / (t1 - t0)
            return tuple(c0[k] + (c1[k] - c0[k]) * s for k in range(3))
    return stops[-1][1]


# . palettes (2026-10-24): the colour maps of util/colormaps.py (the ones of the PES Analysis
#   window, sampled the same way as util/easyplot), plus the built-in default. Diverging maps
#   are used for difference plots (centred on zero); with a sequential choice those keep the
#   built-in blue - white - red.
DEFAULT_PALETTE = "default"
DIVERGING = {"bam", "berlin", "vik", "coolwarm", "coolwarm2"}


def _color_maps():
    try:
        from util.colormaps import COLOR_MAPS
    except ImportError:
        return {}
    return {name.strip(): stops for name, stops in COLOR_MAPS.items()}


_MAPS = _color_maps()


def palette_names():
    """ [ (id, label) ] for a combo box; the first is the built-in default. """
    return [(DEFAULT_PALETTE, "Default")] + [(n, n) for n in _MAPS]


def _sample(stops, t):
    """ As util/easyplot/color_utils.get_color: linear between the stops of the dict. """
    keys = list(stops.keys())
    vals = list(stops.values())
    if t <= keys[0]:
        return tuple(vals[0])
    for k in range(1, len(keys)):
        if t <= keys[k]:
            s = (t - keys[k - 1]) / (keys[k] - keys[k - 1])
            return tuple(vals[k - 1][c] + s * (vals[k][c] - vals[k - 1][c]) for c in range(3))
    return tuple(vals[-1])


def palette_color(t, name=DEFAULT_PALETTE, diverging=False, reverse=False):
    """ t in [0, 1] -> (r, g, b) with the named palette. """
    t = min(1.0, max(0.0, t))
    if reverse:
        t = 1.0 - t
    stops = _MAPS.get(name)
    if stops is None or (diverging and name not in DIVERGING):
        return colormap(t, diverging)
    return _sample(stops, t)


class AnalysisPlot(Gtk.DrawingArea):

    def __init__(self):
        super().__init__()
        self.kind = None
        self.title = ""
        self.xlabel = ""
        self.ylabel = ""
        self.hover = None
        self.on_click = None
        self.data = {}
        self.palette = DEFAULT_PALETTE
        self.palette_reverse = False
        self.add_events(Gdk.EventMask.POINTER_MOTION_MASK | Gdk.EventMask.LEAVE_NOTIFY_MASK |
                        Gdk.EventMask.BUTTON_PRESS_MASK)
        self.connect("draw", self._on_draw)
        self.connect("motion-notify-event", self._on_motion)
        self.connect("leave-notify-event", lambda *_: self._set_hover(None))
        self.connect("button-press-event", self._on_button)
        self.set_size_request(480, 300)

    def set_palette(self, name=DEFAULT_PALETTE, reverse=False):
        """ Colour map of heatmaps, colour bars and coloured scatter points (see palette_names()). """
        self.palette = name or DEFAULT_PALETTE
        self.palette_reverse = bool(reverse)
        self.queue_draw()

    def _color(self, t, diverging=False):
        return palette_color(t, self.palette, diverging, self.palette_reverse)

    # ------------------------------------------------------------------ data
    def clear(self):
        self.kind = None
        self.data = {}
        self.hover = None
        self.queue_draw()

    def set_lines(self, series, xlabel="", ylabel="", title="", ymin=None, hline=None, xlabels=None):
        for k, s in enumerate(series):
            s.setdefault("color", PALETTE[k % len(PALETTE)])
            s.setdefault("style", "line")
            s["x"] = [float(v) for v in s["x"]]
            s["y"] = [float(v) for v in s["y"]]
        self.kind, self.title, self.xlabel, self.ylabel = "lines", title, xlabel, ylabel
        self.data = {"series": series, "ymin": ymin, "hline": hline, "xlabels": xlabels}
        self.hover = None
        self.queue_draw()

    def set_heatmap(self, matrix, xlabel="", ylabel="", title="", vmin=None, vmax=None, diverging=False,
                    row_labels=None, col_labels=None, value_label="value", overlays=None):
        """overlays (2026-10-23): [{"cells": [(row, col), ...] (cell units, floats allowed),
        "style": "box" | "path" | "dots", "color": (r, g, b), "width": px}] drawn on top.
        nan cells are grey (no data)."""
        rows = len(matrix)
        cols = len(matrix[0]) if rows else 0
        flat = [v for row in matrix for v in row if v == v]
        lo = (min(flat) if flat else 0.0) if vmin is None else vmin
        hi = (max(flat) if flat else 1.0) if vmax is None else vmax
        if hi == lo:
            hi = lo + 1.0
        self.kind, self.title, self.xlabel, self.ylabel = "heatmap", title, xlabel, ylabel
        self.data = {"m": matrix, "rows": rows, "cols": cols, "lo": lo, "hi": hi, "div": diverging,
                     "rl": row_labels, "cl": col_labels, "vlabel": value_label, "overlays": overlays or []}
        self.hover = None
        self.queue_draw()

    def set_categories(self, codes, colors, names, xlabel="", ylabel="", title="", row_labels=None, col_labels=None):
        """codes[row][col] -> one key of colors. Rows are drawn top to bottom."""
        self.kind, self.title, self.xlabel, self.ylabel = "categories", title, xlabel, ylabel
        self.data = {"m": codes, "rows": len(codes), "cols": len(codes[0]) if len(codes) else 0,
                     "colors": colors, "names": names, "rl": row_labels, "cl": col_labels}
        self.hover = None
        self.queue_draw()

    def set_scatter(self, x, y, c=None, xlabel="", ylabel="", title="", clabel=""):
        self.kind, self.title, self.xlabel, self.ylabel = "scatter", title, xlabel, ylabel
        self.data = {"x": [float(v) for v in x], "y": [float(v) for v in y],
                     "c": None if c is None else [float(v) for v in c], "clabel": clabel}
        self.hover = None
        self.queue_draw()

    # ------------------------------------------------------------------ geometry
    def _frame(self, width, height):
        right = 90 if self.kind in ("heatmap", "scatter") and (self.kind != "scatter" or self.data.get("c")) else 20
        if self.kind == "categories":
            right = 120
        if self.kind == "lines" and len(self.data.get("series", [])) > 1:
            right = 20
        top = 34 if self.title else 16
        return 72, top, max(width - 72 - right, 20), max(height - top - 56, 20)

    def _line_ranges(self):
        xs = [v for s in self.data["series"] for v in s["x"]]
        ys = [v for s in self.data["series"] for v in s["y"] if v == v]
        if not xs or not ys:
            return 0, 1, 0, 1
        x0, x1 = min(xs), max(xs)
        y0, y1 = min(ys), max(ys)
        if any(s["style"] == "bars" for s in self.data["series"]):
            y0 = min(0.0, y0)
            x0 -= 0.5; x1 += 0.5
        if self.data.get("ymin") is not None:
            y0 = min(self.data["ymin"], y0)
        if x1 - x0 <= 1e-9 * max(1.0, abs(x1)):
            x0, x1 = x0 - 1, x1 + 1
        # . 2026-10-24: a (nearly) constant series -- e.g. the density of a box of fixed volume,
        #   equal to ~1e-16 -- gave a tick step below the float resolution: endless tick loop
        if y1 - y0 <= 1e-9 * max(1.0, abs(y1)):
            pad = 0.06 * (abs(y1) or 1.0)
            return x0, x1, y0 - pad, y1 + pad
        pad = 0.06 * (y1 - y0 or abs(y1) or 1.0)
        return x0, x1, y0 - (pad if self.data.get("ymin") is None else 0), y1 + pad

    # ------------------------------------------------------------------ drawing
    def _on_draw(self, widget, cr, width=None, height=None):
        width = width or widget.get_allocated_width()
        height = height or widget.get_allocated_height()
        cr.set_source_rgb(1, 1, 1)
        cr.paint()
        if self.kind is None:
            return False
        cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        x0, y0, w, h = self._frame(width, height)
        getattr(self, "_draw_" + self.kind)(cr, x0, y0, w, h)
        # . titles
        cr.set_source_rgb(0.1, 0.1, 0.1)
        cr.set_font_size(12)
        if self.xlabel:
            e = cr.text_extents(self.xlabel)
            cr.move_to(x0 + w / 2 - e.width / 2, y0 + h + 42); cr.show_text(self.xlabel)
        if self.ylabel:
            cr.save(); e = cr.text_extents(self.ylabel)
            cr.move_to(16, y0 + h / 2 + e.width / 2); cr.rotate(-math.pi / 2); cr.show_text(self.ylabel); cr.restore()
        if self.title:
            cr.set_font_size(13)
            e = cr.text_extents(self.title)
            cr.move_to(x0 + w / 2 - e.width / 2, y0 - 12); cr.show_text(self.title)
        if self.hover is not None and self.hover.get("text"):
            text = self.hover["text"]
            cr.set_font_size(11)
            e = cr.text_extents(text)
            bx, by = max(x0 + w - e.width - 14, 4), y0 + h - e.height - 14
            cr.set_source_rgba(1, 1, 1, 0.92); cr.rectangle(bx - 5, by - 2, e.width + 12, e.height + 9); cr.fill()
            cr.set_source_rgb(0.3, 0.3, 0.3); cr.set_line_width(0.6); cr.rectangle(bx - 5, by - 2, e.width + 12, e.height + 9); cr.stroke()
            cr.set_source_rgb(0.1, 0.1, 0.1); cr.move_to(bx, by + e.height + 2); cr.show_text(text)
        return False

    def _axes(self, cr, x0, y0, w, h, xr, yr, xlabels=None):
        cr.set_font_size(10.5)
        cr.set_line_width(1)
        (a0, a1), (b0, b1) = xr, yr
        if xlabels is None:
            for v in _ticks(a0, a1, 8):
                x = x0 + (v - a0) / (a1 - a0) * w
                cr.set_source_rgb(0.92, 0.92, 0.92); cr.move_to(x, y0); cr.line_to(x, y0 + h); cr.stroke()
                cr.set_source_rgb(0.2, 0.2, 0.2)
                t = _fmt(v); e = cr.text_extents(t); cr.move_to(x - e.width / 2, y0 + h + 15); cr.show_text(t)
        else:
            every = max(1, int(len(xlabels) * 7 / max(w, 1)) + 1)
            for k, label in enumerate(xlabels):
                if k % every:
                    continue
                x = x0 + (k - a0) / (a1 - a0) * w
                cr.save(); cr.set_source_rgb(0.2, 0.2, 0.2)
                cr.move_to(x + 3, y0 + h + 6); cr.rotate(math.pi / 4); cr.show_text(str(label)); cr.restore()
        for v in _ticks(b0, b1, 6):
            y = y0 + h - (v - b0) / (b1 - b0) * h
            cr.set_source_rgb(0.92, 0.92, 0.92); cr.move_to(x0, y); cr.line_to(x0 + w, y); cr.stroke()
            cr.set_source_rgb(0.2, 0.2, 0.2)
            t = _fmt(v); e = cr.text_extents(t); cr.move_to(x0 - e.width - 6, y + e.height / 2); cr.show_text(t)

    def _draw_lines(self, cr, x0, y0, w, h):
        a0, a1, b0, b1 = self._line_ranges()
        self._axes(cr, x0, y0, w, h, (a0, a1), (b0, b1), self.data.get("xlabels"))
        X = lambda v: x0 + (v - a0) / (a1 - a0) * w
        Y = lambda v: y0 + h - (v - b0) / (b1 - b0) * h
        cr.save(); cr.rectangle(x0, y0, w, h); cr.clip()
        if self.data.get("hline") is not None:
            cr.set_source_rgb(0.55, 0.55, 0.55); cr.set_dash([4, 4]); cr.set_line_width(1)
            cr.move_to(x0, Y(self.data["hline"])); cr.line_to(x0 + w, Y(self.data["hline"])); cr.stroke(); cr.set_dash([])
        series = self.data["series"]
        nbars = sum(1 for s in series if s["style"] == "bars")
        bar_k = 0
        for s in series:
            cr.set_source_rgb(*s["color"])
            if s["style"] == "bars":
                bw = w / max(a1 - a0, 1) * 0.8 / max(nbars, 1)
                colors = s.get("colors")             # . optional colour per bar
                for m, (x, y) in enumerate(zip(s["x"], s["y"])):
                    if y != y:
                        continue
                    bx = X(x) - 0.4 * w / max(a1 - a0, 1) + bar_k * bw
                    if colors:
                        cr.set_source_rgb(*colors[m])
                    cr.rectangle(bx, min(Y(y), Y(0)), max(bw, 1), abs(Y(y) - Y(0)))
                    if colors:
                        cr.fill()
                cr.fill()
                bar_k += 1
            elif s["style"] == "points":
                for x, y in zip(s["x"], s["y"]):
                    cr.arc(X(x), Y(y), 2.2, 0, 2 * math.pi); cr.fill()
            else:
                cr.set_line_width(s.get("width", 1.6))
                if s.get("dash"):
                    cr.set_dash([6, 4])
                first = True
                for x, y in zip(s["x"], s["y"]):
                    if y != y:
                        first = True
                        continue
                    (cr.move_to if first else cr.line_to)(X(x), Y(y))
                    first = False
                cr.stroke(); cr.set_dash([])
        if self.hover is not None and "pt" in self.hover:
            px, py = self.hover["pt"]
            cr.set_source_rgb(0.85, 0.1, 0.1); cr.arc(X(px), Y(py), 3.6, 0, 2 * math.pi); cr.fill()
        cr.restore()
        cr.set_source_rgb(0.2, 0.2, 0.2); cr.set_line_width(1.1); cr.rectangle(x0, y0, w, h); cr.stroke()
        # . legend
        labelled = [s for s in series if s.get("label")]
        if len(labelled) > 1:
            cr.set_font_size(10.5)
            lw = max(cr.text_extents(s["label"]).width for s in labelled) + 34
            lx, ly = x0 + w - lw - 6, y0 + 10
            cr.set_source_rgba(1, 1, 1, 0.88)
            cr.rectangle(lx - 6, ly - 6, lw, 15 * len(labelled) + 6); cr.fill()
            for s in labelled:
                cr.set_source_rgb(*s["color"]); cr.rectangle(lx, ly, 14, 4); cr.fill()
                cr.set_source_rgb(0.15, 0.15, 0.15); cr.move_to(lx + 20, ly + 6); cr.show_text(s["label"])
                ly += 15

    def _colorbar(self, cr, x, y0, h, lo, hi, diverging, label):
        steps = 60
        for k in range(steps):
            t = 1 - k / steps
            cr.set_source_rgb(*self._color(t, diverging))
            cr.rectangle(x, y0 + k * h / steps, 14, h / steps + 1); cr.fill()
        cr.set_source_rgb(0.2, 0.2, 0.2); cr.set_line_width(0.8); cr.rectangle(x, y0, 14, h); cr.stroke()
        cr.set_font_size(10)
        for v in _ticks(lo, hi, 5):
            yy = y0 + h - (v - lo) / (hi - lo) * h
            cr.move_to(x + 18, yy + 4); cr.show_text(_fmt(v))
        if label:
            cr.save(); cr.move_to(x + 70, y0 + h / 2 + cr.text_extents(label).width / 2)
            cr.rotate(-math.pi / 2); cr.show_text(label); cr.restore()

    def _draw_matrix_axes(self, cr, x0, y0, w, h, rows, cols, rl, cl):
        cr.set_font_size(10)
        cr.set_source_rgb(0.2, 0.2, 0.2)
        every_c = max(1, int(cols * (24 if cl else 40) / max(w, 1)) + 1)
        for k in range(0, cols, every_c):
            label = cl[k] if cl else str(k)
            x = x0 + (k + 0.5) / cols * w
            if cl:
                cr.save(); cr.move_to(x + 2, y0 + h + 6); cr.rotate(math.pi / 4); cr.show_text(str(label)); cr.restore()
            else:
                e = cr.text_extents(str(label)); cr.move_to(x - e.width / 2, y0 + h + 14); cr.show_text(str(label))
        every_r = max(1, int(rows * 13 / max(h, 1)) + 1)
        for k in range(0, rows, every_r):
            label = rl[k] if rl else str(k)
            y = y0 + (k + 0.5) / rows * h
            e = cr.text_extents(str(label)); cr.move_to(x0 - e.width - 5, y + 4); cr.show_text(str(label))

    def _draw_heatmap(self, cr, x0, y0, w, h):
        d = self.data
        rows, cols, lo, hi = d["rows"], d["cols"], d["lo"], d["hi"]
        cw, ch = w / cols, h / rows
        for i in range(rows):
            row = d["m"][i]
            for j in range(cols):
                v = row[j]
                if v != v:
                    cr.set_source_rgb(0.86, 0.86, 0.86)
                else:
                    cr.set_source_rgb(*self._color((v - lo) / (hi - lo), d["div"]))
                cr.rectangle(x0 + j * cw, y0 + i * ch, cw + 0.6, ch + 0.6)
                cr.fill()
        for ov in d.get("overlays") or []:
            cells = ov.get("cells") or []
            cr.set_source_rgb(*ov.get("color", (0.85, 0.1, 0.1)))
            cr.set_line_width(ov.get("width", 2.0))
            style = ov.get("style", "box")
            if style == "box":
                for r, c in cells:
                    cr.rectangle(x0 + c * cw + 1, y0 + r * ch + 1, max(cw - 2, 1), max(ch - 2, 1)); cr.stroke()
            elif style == "path" and len(cells) > 1:
                for n, (r, c) in enumerate(cells):
                    (cr.move_to if n == 0 else cr.line_to)(x0 + (c + 0.5) * cw, y0 + (r + 0.5) * ch)
                cr.stroke()
            else:
                rad = ov.get("radius", 3.0)
                for r, c in cells:
                    cr.arc(x0 + (c + 0.5) * cw, y0 + (r + 0.5) * ch, rad, 0, 2 * math.pi); cr.fill()
            if ov.get("label_cells"):
                cr.set_font_size(11)
                for (r, c), text in ov["label_cells"]:
                    cr.move_to(x0 + (c + 0.5) * cw + 4, y0 + (r + 0.5) * ch - 4); cr.show_text(text)
        cr.set_source_rgb(0.2, 0.2, 0.2); cr.set_line_width(1.1); cr.rectangle(x0, y0, w, h); cr.stroke()
        self._draw_matrix_axes(cr, x0, y0, w, h, rows, cols, d["rl"], d["cl"])
        self._colorbar(cr, x0 + w + 12, y0, h, lo, hi, d["div"], d["vlabel"])
        if self.hover is not None and "cell" in self.hover:
            i, j = self.hover["cell"]
            cr.set_source_rgb(0, 0, 0); cr.set_line_width(1.5)
            cr.rectangle(x0 + j * cw, y0 + i * ch, max(cw, 2), max(ch, 2)); cr.stroke()

    def _draw_categories(self, cr, x0, y0, w, h):
        d = self.data
        rows, cols = d["rows"], d["cols"]
        cw, ch = w / cols, h / rows
        for i in range(rows):
            row = d["m"][i]
            for j in range(cols):
                cr.set_source_rgb(*d["colors"].get(row[j], (1, 1, 1)))
                cr.rectangle(x0 + j * cw, y0 + i * ch, cw + 0.6, ch + 0.6)
                cr.fill()
        cr.set_source_rgb(0.2, 0.2, 0.2); cr.set_line_width(1.1); cr.rectangle(x0, y0, w, h); cr.stroke()
        self._draw_matrix_axes(cr, x0, y0, w, h, rows, cols, d["rl"], d["cl"])
        cr.set_font_size(10.5)
        ly = y0
        for key, name in d["names"].items():
            cr.set_source_rgb(*d["colors"].get(key, (1, 1, 1))); cr.rectangle(x0 + w + 12, ly, 14, 12); cr.fill()
            cr.set_source_rgb(0.3, 0.3, 0.3); cr.set_line_width(0.6); cr.rectangle(x0 + w + 12, ly, 14, 12); cr.stroke()
            cr.set_source_rgb(0.15, 0.15, 0.15); cr.move_to(x0 + w + 32, ly + 10); cr.show_text(name)
            ly += 17

    def _draw_scatter(self, cr, x0, y0, w, h):
        d = self.data
        if not d["x"]:
            return
        a0, a1 = min(d["x"]), max(d["x"]); b0, b1 = min(d["y"]), max(d["y"])
        pa, pb = 0.06 * (a1 - a0 or 1), 0.06 * (b1 - b0 or 1)
        a0, a1, b0, b1 = a0 - pa, a1 + pa, b0 - pb, b1 + pb
        self._axes(cr, x0, y0, w, h, (a0, a1), (b0, b1))
        X = lambda v: x0 + (v - a0) / (a1 - a0) * w
        Y = lambda v: y0 + h - (v - b0) / (b1 - b0) * h
        c = d["c"]
        if c:
            c0, c1 = min(c), max(c)
            c1 = c1 if c1 > c0 else c0 + 1
        for k, (x, y) in enumerate(zip(d["x"], d["y"])):
            cr.set_source_rgb(*(self._color((c[k] - c0) / (c1 - c0)) if c else PALETTE[0]))
            cr.arc(X(x), Y(y), 3.0, 0, 2 * math.pi); cr.fill()
        if self.hover is not None and "pt" in self.hover:
            px, py = self.hover["pt"]
            cr.set_source_rgb(0.85, 0.1, 0.1); cr.set_line_width(1.6); cr.arc(X(px), Y(py), 5, 0, 2 * math.pi); cr.stroke()
        cr.set_source_rgb(0.2, 0.2, 0.2); cr.set_line_width(1.1); cr.rectangle(x0, y0, w, h); cr.stroke()
        if c:
            self._colorbar(cr, x0 + w + 12, y0, h, c0, c1, False, d["clabel"])

    # ------------------------------------------------------------------ interaction
    def _locate(self, ex, ey):
        width, height = self.get_allocated_width(), self.get_allocated_height()
        x0, y0, w, h = self._frame(width, height)
        if not (x0 <= ex <= x0 + w and y0 <= ey <= y0 + h):
            return None
        d = self.data
        if self.kind == "lines":
            a0, a1, b0, b1 = self._line_ranges()
            best = None
            for s in d["series"]:
                for k, (x, y) in enumerate(zip(s["x"], s["y"])):
                    if y != y:
                        continue
                    px = x0 + (x - a0) / (a1 - a0) * w
                    py = y0 + h - (y - b0) / (b1 - b0) * h
                    dist = abs(px - ex) + (0.25 * abs(py - ey) if s["style"] != "bars" else 0)
                    if best is None or dist < best[0]:
                        best = (dist, s, k, x, y)
            if best is None:
                return None
            _, s, k, x, y = best
            xl = d.get("xlabels")
            xt = xl[int(round(x))] if xl and 0 <= int(round(x)) < len(xl) else _fmt(x)
            text = "{}{} = {},  {} = {}".format((s.get("label") + ":  ") if s.get("label") else "",
                                                 self.xlabel.split(" (")[0] or "x", xt,
                                                 self.ylabel.split(" (")[0] or "y", _fmt(y))
            return {"pt": (x, y), "text": text, "x": x, "y": y, "series": s.get("label"), "index": k}
        if self.kind in ("heatmap", "categories"):
            j = int((ex - x0) / w * d["cols"]); i = int((ey - y0) / h * d["rows"])
            i, j = min(i, d["rows"] - 1), min(j, d["cols"] - 1)
            rl = d["rl"][i] if d["rl"] else i
            cl = d["cl"][j] if d["cl"] else j
            v = d["m"][i][j]
            vt = d["names"].get(v, v) if self.kind == "categories" else _fmt(v)
            return {"cell": (i, j), "text": "{} / {}:  {}".format(rl, cl, vt), "row": i, "col": j}
        if self.kind == "scatter":
            a0, a1 = min(d["x"]), max(d["x"]); b0, b1 = min(d["y"]), max(d["y"])
            pa, pb = 0.06 * (a1 - a0 or 1), 0.06 * (b1 - b0 or 1)
            a0, a1, b0, b1 = a0 - pa, a1 + pa, b0 - pb, b1 + pb
            best = None
            for k, (x, y) in enumerate(zip(d["x"], d["y"])):
                px = x0 + (x - a0) / (a1 - a0) * w
                py = y0 + h - (y - b0) / (b1 - b0) * h
                dist = (px - ex) ** 2 + (py - ey) ** 2
                if best is None or dist < best[0]:
                    best = (dist, k, x, y)
            if best is None or best[0] > 400:
                return None
            _, k, x, y = best
            extra = "" if not d["c"] else ",  {} = {}".format(d["clabel"] or "c", _fmt(d["c"][k]))
            return {"pt": (x, y), "text": "{} = {},  {} = {}{}".format(self.xlabel or "x", _fmt(x), self.ylabel or "y",
                                                                     _fmt(y), extra), "index": k}
        return None

    def _set_hover(self, hover):
        if hover != self.hover:
            self.hover = hover
            self.queue_draw()
        return False

    def _on_motion(self, _widget, event):
        if self.kind is None:
            return False
        return self._set_hover(self._locate(event.x, event.y))

    def _on_button(self, _widget, event):
        if event.button != 1 or self.on_click is None or self.kind is None:
            return False
        info = self._locate(event.x, event.y)
        if info is not None:
            self.on_click(info)
        return True

    def export(self, path, width=1000, height=620):
        if path.lower().endswith(".svg"):
            surface = cairo.SVGSurface(path, width, height)
        else:
            surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        cr = cairo.Context(surface)
        hover, self.hover = self.hover, None
        self._on_draw(self, cr, width, height)
        self.hover = hover
        if path.lower().endswith(".svg"):
            surface.finish()
        else:
            surface.write_to_png(path)


def show_plot_window(window_title, series, xlabel="", ylabel="", title="", ymin=None, default_name="plot.png",
                     parent=None, size=(760, 470)):
    """ 2026-10-24: a window with one AnalysisPlot of line series (axis titles, hover readout)
        and an "Export Image..." button (PNG or SVG). Used by the RMSD, RMSF and Box Density
        windows (they drew with the old XYPlot: no axis titles, negative y axis, no export).
        series: as AnalysisPlot.set_lines(). Returns ( window, plot ). """
    window = Gtk.Window(title=window_title)
    if parent is not None:
        window.set_transient_for(parent)
    box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4, margin=6)
    plot = AnalysisPlot()
    plot.set_lines(series, xlabel, ylabel, title, ymin=ymin)
    box.pack_start(plot, True, True, 0)
    button = Gtk.Button(label="Export Image...")
    button.set_halign(Gtk.Align.END)

    def on_export(_b):
        dialog = Gtk.FileChooserDialog(title="Export Image", parent=window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name(default_name)
        response, path = dialog.run(), dialog.get_filename()
        dialog.destroy()
        if response == Gtk.ResponseType.OK and path:
            if not path.lower().endswith((".png", ".svg")):
                path += ".png"
            plot.export(path)
    button.connect("clicked", on_export)
    box.pack_start(button, False, False, 0)
    window.add(box)
    window.set_default_size(*size)
    window.show_all()
    return window, plot
