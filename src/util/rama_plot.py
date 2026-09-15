#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Ramachandran plot widget
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      RamachandranPlot: a plain Cairo GtkDrawingArea (same family as
#      util/easyplot's XYPlot/ImagePlot -- on_draw(widget, cr,
#      export_scale=1), so util.easyplot.export_utils.export_plot_to_png()
#      works on it unmodified) specialised for phi/psi scatter plots.
#      Kept as its own widget (rather than reusing XYPlot/XYScatterPlot)
#      because a Ramachandran plot has requirements neither of those
#      generic classes cover: a fixed -180..180 SQUARE domain (never
#      auto-scaled to the data), shaded background regions, and points
#      coloured by residue category rather than by series.
#
#      Fed by util/ramachandran.py's compute_phi_psi() (GTK-free); this
#      module only draws, it does not touch pDynamo/vismol at all.
#
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
import cairo
import math

# . Category -> RGB fill colour. Kept as a module-level dict (rather than
#   an attribute) so ramachandran_analysis_window.py can reuse the exact
#   same colours when building its own legend/status text.
CATEGORY_COLORS = {
    "general":     (0.20, 0.40, 0.80),
    "glycine":     (0.95, 0.60, 0.15),
    "proline":     (0.25, 0.65, 0.35),
    "pre-proline": (0.55, 0.35, 0.75),
}
CATEGORY_LABELS = {
    "general":     "Geral",
    "glycine":     "Glicina",
    "proline":     "Prolina",
    "pre-proline": "Pré-prolina",
}
# . Order in which categories are drawn (and listed in the legend) --
#   "general" first/bottom so the smaller, more informative categories
#   (Gly/Pro/pre-Pro) are never hidden underneath a sea of generic points.
_CATEGORY_ORDER = ["general", "glycine", "proline", "pre-proline"]


def _ellipse_points(cx, cy, rx, ry, n=48):
    """ n points around an ellipse centred at (cx, cy) in (phi, psi)
        degrees -- used to approximate the classic Ramachandran favoured
        regions with smooth shapes instead of hand-typed polygons.
    """
    points = []
    for i in range(n):
        t = 2.0 * math.pi * i / n
        points.append((cx + rx * math.cos(t), cy + ry * math.sin(t)))
    return points


# . Illustrative approximations of the classic Ramachandran regions,
#   each with a smaller "core favoured" ellipse and a larger, lighter
#   "generously allowed" ellipse around it. These are NOT statistical
#   contours fitted to real structure databases (unlike e.g. the
#   MolProbity/Top500 contours) -- they exist to give a qualitative
#   sense of where a point falls at a glance. Do not use this plot as
#   the sole basis for a structure-validation judgement; for that,
#   cross-check against a dedicated validation tool.
#
#   Split by category rather than one-size-fits-all, because the two
#   most common non-generic residues have genuinely different backbone
#   geometry, not just different plotting conventions:
#
#   - Glycine has no side chain, so its two Halpha hydrogens are
#     equivalent -- its Ramachandran map is symmetric under
#     (phi, psi) -> (-phi, -psi). _GLYCINE_EXTRA_REGIONS are exactly the
#     point-mirror images of the general alpha-R and beta regions, drawn
#     ON TOP OF (in addition to, never instead of) the general regions
#     below -- a glycine is still free to sit in a "general" region too.
#   - Proline's pyrrolidine ring covalently closes back onto the
#     backbone nitrogen, locking phi to approximately -60 degrees
#     regardless of psi -- so _PROLINE_REGIONS is a narrow vertical
#     band (tight in phi, generous in psi) instead of the general
#     regions, which is why it is drawn INSTEAD OF them (see on_draw()):
#     a proline cannot reach the general alpha-L/beta-mirror areas that
#     glycine can, so showing the general background for it would
#     overstate what is actually allowed.
_GENERAL_REGIONS = [
    {"label": "αR", "core": (-63, -43, 26, 27), "allowed": (-63, -43, 48, 48)},
    {"label": "β",  "core": (-120, 130, 55, 40), "allowed": (-120, 130, 85, 60)},
    {"label": "L",       "core": (57, 40, 20, 20),    "allowed": (57, 40, 35, 35)},
]

_GLYCINE_EXTRA_REGIONS = [
    # Mirror of alpha-R: near-coincides with the general "L" spot above,
    # just reinforcing it -- no separate label, "L" already marks it.
    {"label": None,  "core": (63, 43, 26, 27),     "allowed": (63, 43, 48, 48)},
    {"label": "Gly", "core": (120, -130, 55, 40),  "allowed": (120, -130, 85, 60)},   # mirror of beta -- Gly-only territory
]

_PROLINE_REGIONS = [
    # Alpha-like: near-coincides with the general alpha-R region above,
    # just narrower in phi -- no separate label, "αR" already marks it.
    {"label": None,  "core": (-60, -30, 16, 26), "allowed": (-60, -30, 26, 36)},
    {"label": "Pro", "core": (-60, 145, 16, 32), "allowed": (-60, 145, 26, 46)},   # polyproline-II/extended
]

_REGION_CORE_COLOR    = (0.55, 0.80, 0.55, 0.55)
_REGION_ALLOWED_COLOR = (0.55, 0.80, 0.55, 0.25)
_REGION_LABEL_COLOR   = (0.20, 0.40, 0.20)

_POINT_RADIUS = 3.0
_HOVER_PICK_RADIUS_PX = 7.0


class RamachandranPlot(Gtk.DrawingArea):
    """ Ramachandran (phi/psi) scatter plot.

        Usage:
            plot = RamachandranPlot()
            plot.set_data(results)   # results: see util.ramachandran.compute_phi_psi()
            plot.on_hover_callback = my_handler   # optional, see below

        on_hover_callback, if set, is called with either one point-dict
        (the nearest point to the mouse, same shape as a compute_phi_psi()
        row, when the cursor is within _HOVER_PICK_RADIUS_PX pixels of a
        plotted point) or None (cursor not near any point / left the
        widget) -- lets the containing window show "residue under the
        mouse" feedback without this widget needing to know anything
        about GTK labels or window layout.
    """

    def __init__(self):
        super().__init__()
        self.data = []
        self._hit_test_cache = []  # [(canvas_x, canvas_y, point_dict), ...], rebuilt on_draw
        self.on_hover_callback = None

        self.add_events(Gdk.EventMask.POINTER_MOTION_MASK |
                         Gdk.EventMask.LEAVE_NOTIFY_MASK)
        self.connect("draw", self.on_draw)
        self.connect("motion-notify-event", self._on_motion)
        self.connect("leave-notify-event", self._on_leave)

    def set_data(self, points):
        """ points: list of dicts as returned by
            util.ramachandran.compute_phi_psi().
        """
        self.data = points
        self.queue_draw()

    def _layout(self, width, height):
        """ Returns (origin_x, origin_y, size) for the plot's square
            drawing area (in "phi -180..180 maps to [origin_x, origin_x
            + size]" pixel space), centred within the widget with room
            left for axis labels.
        """
        left, right, top, bottom = 56, 22, 30, 46
        avail_w = max(width  - left - right, 10)
        avail_h = max(height - top  - bottom, 10)
        size = min(avail_w, avail_h)
        origin_x = left   + (avail_w - size) / 2.0
        origin_y = top    + (avail_h - size) / 2.0
        return origin_x, origin_y, size

    def _to_canvas(self, phi, psi, origin_x, origin_y, size):
        x = origin_x + ((phi + 180.0) / 360.0) * size
        y = origin_y + size - ((psi + 180.0) / 360.0) * size
        return x, y

    def _from_canvas(self, x, y, origin_x, origin_y, size):
        phi = (x - origin_x) / size * 360.0 - 180.0
        psi = (origin_y + size - y) / size * 360.0 - 180.0
        return phi, psi

    def _draw_region(self, cr, origin_x, origin_y, size, ellipse, color):
        cx, cy, rx, ry = ellipse
        cr.set_source_rgba(*color)
        first = True
        for phi, psi in _ellipse_points(cx, cy, rx, ry):
            x, y = self._to_canvas(phi, psi, origin_x, origin_y, size)
            if first:
                cr.move_to(x, y)
                first = False
            else:
                cr.line_to(x, y)
        cr.close_path()
        cr.fill()

    def _draw_region_label(self, cr, origin_x, origin_y, size, region):
        label = region.get("label")
        if not label:
            return
        cx, cy, _rx, _ry = region["core"]
        x, y = self._to_canvas(cx, cy, origin_x, origin_y, size)
        cr.select_font_face("Sans", cairo.FONT_SLANT_ITALIC, cairo.FONT_WEIGHT_NORMAL)
        cr.set_font_size(11)
        cr.set_source_rgb(*_REGION_LABEL_COLOR)
        extents = cr.text_extents(label)
        cr.move_to(x - extents.width / 2.0, y + extents.height / 2.0)
        cr.show_text(label)

    def on_draw(self, widget, cr, export_scale=1):
        width  = widget.get_allocated_width()
        height = widget.get_allocated_height()
        if export_scale != 1:
            cr.scale(export_scale, export_scale)

        cr.set_source_rgb(1, 1, 1)
        cr.paint()

        origin_x, origin_y, size = self._layout(width, height)
        step = size / 6.0  # 60 degrees per major step

        by_category = {}
        for point in self.data:
            by_category.setdefault(point.get("category", "general"), []).append(point)

        # . The general regions stay on the plot regardless (the general
        #   category is almost always the majority of a real protein's
        #   residues); glycine's and proline's regions are ADDED on top
        #   of that, only once that category actually has a point in the
        #   current data -- so an all-general system is not cluttered
        #   with a Pro/Gly band that has nothing to do with it.
        regions_to_draw = list(_GENERAL_REGIONS)
        if "glycine" in by_category:
            regions_to_draw += _GLYCINE_EXTRA_REGIONS
        if "proline" in by_category:
            regions_to_draw += _PROLINE_REGIONS

        cr.save()
        cr.rectangle(origin_x, origin_y, size, size)
        cr.clip()

        # . Minor grid every 30 degrees, light grey.
        cr.set_line_width(1)
        cr.set_source_rgb(0.88, 0.88, 0.88)
        for i in range(13):
            x = origin_x + i * (size / 12.0)
            y = origin_y + i * (size / 12.0)
            cr.move_to(x, origin_y); cr.line_to(x, origin_y + size)
            cr.move_to(origin_x, y); cr.line_to(origin_x + size, y)
        cr.stroke()

        # . Favoured regions (approximate, see _GENERAL_REGIONS docstring above).
        for region in regions_to_draw:
            self._draw_region(cr, origin_x, origin_y, size, region["allowed"], _REGION_ALLOWED_COLOR)
        for region in regions_to_draw:
            self._draw_region(cr, origin_x, origin_y, size, region["core"], _REGION_CORE_COLOR)
        for region in regions_to_draw:
            self._draw_region_label(cr, origin_x, origin_y, size, region)

        # . Major grid every 60 degrees, darker.
        cr.set_line_width(1)
        cr.set_source_rgb(0.7, 0.7, 0.7)
        for i in range(7):
            x = origin_x + i * step
            y = origin_y + i * step
            cr.move_to(x, origin_y); cr.line_to(x, origin_y + size)
            cr.move_to(origin_x, y); cr.line_to(origin_x + size, y)
        cr.stroke()

        # . phi = 0 / psi = 0 axes.
        cr.set_source_rgb(0.15, 0.15, 0.15)
        cr.set_line_width(1.6)
        cx, cy = self._to_canvas(0, 0, origin_x, origin_y, size)
        cr.move_to(cx, origin_y); cr.line_to(cx, origin_y + size)
        cr.move_to(origin_x, cy); cr.line_to(origin_x + size, cy)
        cr.stroke()

        # . Scatter points, drawn category-by-category so the order in
        #   _CATEGORY_ORDER controls what ends up on top.
        self._hit_test_cache = []
        for category in _CATEGORY_ORDER:
            color = CATEGORY_COLORS.get(category, CATEGORY_COLORS["general"])
            cr.set_source_rgb(*color)
            for point in by_category.get(category, []):
                x, y = self._to_canvas(point["phi"], point["psi"], origin_x, origin_y, size)
                cr.arc(x, y, _POINT_RADIUS, 0, 2 * math.pi)
                cr.fill()
                self._hit_test_cache.append((x, y, point))

        cr.restore()  # end clip to the plot square

        # . Frame around the square.
        cr.set_source_rgb(0.15, 0.15, 0.15)
        cr.set_line_width(1.6)
        cr.rectangle(origin_x, origin_y, size, size)
        cr.stroke()

        # . Tick labels every 60 degrees.
        cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        cr.set_font_size(10)
        cr.set_source_rgb(0.2, 0.2, 0.2)
        for i in range(7):
            label = str(-180 + i * 60)
            x = origin_x + i * step
            extents = cr.text_extents(label)
            cr.move_to(x - extents.width / 2.0, origin_y + size + 16)
            cr.show_text(label)
            y = origin_y + size - i * step
            cr.move_to(origin_x - 30, y + 3)
            cr.show_text(label)

        # . Axis titles.
        cr.set_font_size(12)
        phi_label = "φ (graus)"
        extents = cr.text_extents(phi_label)
        cr.move_to(origin_x + size / 2.0 - extents.width / 2.0, origin_y + size + 34)
        cr.show_text(phi_label)

        cr.save()
        cr.translate(16, origin_y + size / 2.0)
        cr.rotate(-math.pi / 2.0)
        psi_label = "ψ (graus)"
        extents = cr.text_extents(psi_label)
        cr.move_to(-extents.width / 2.0, 0)
        cr.show_text(psi_label)
        cr.restore()

        self._draw_legend(cr, origin_x, origin_y, size, by_category)
        return False

    def _draw_legend(self, cr, origin_x, origin_y, size, by_category):
        present = [c for c in _CATEGORY_ORDER if by_category.get(c)]
        if not present:
            return

        # . Bottom-right corner of the square (large positive phi, very
        #   negative psi) is essentially always empty for real proteins
        #   -- the same reason that corner is coloured "disallowed" on
        #   any textbook Ramachandran plot -- so the legend box never
        #   covers a favoured region or real data.
        row_h = 15
        box_w = 118
        box_h = 10 + row_h * len(present)
        box_x = origin_x + size - box_w - 8
        box_y = origin_y + size - box_h - 8

        cr.set_source_rgba(1, 1, 1, 0.85)
        cr.rectangle(box_x, box_y, box_w, box_h)
        cr.fill()
        cr.set_source_rgba(0.6, 0.6, 0.6, 0.8)
        cr.set_line_width(1)
        cr.rectangle(box_x, box_y, box_w, box_h)
        cr.stroke()

        cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        cr.set_font_size(10)
        for i, category in enumerate(present):
            cy = box_y + 10 + i * row_h + 5
            cr.set_source_rgb(*CATEGORY_COLORS[category])
            cr.arc(box_x + 12, cy, 4, 0, 2 * math.pi)
            cr.fill()
            cr.set_source_rgb(0.15, 0.15, 0.15)
            cr.move_to(box_x + 22, cy + 3.5)
            cr.show_text(CATEGORY_LABELS[category])

    def _on_motion(self, widget, event):
        if not self._hit_test_cache:
            return False
        nearest = None
        nearest_dist2 = (_HOVER_PICK_RADIUS_PX * _HOVER_PICK_RADIUS_PX)
        for x, y, point in self._hit_test_cache:
            dx = x - event.x
            dy = y - event.y
            dist2 = dx * dx + dy * dy
            if dist2 <= nearest_dist2:
                nearest_dist2 = dist2
                nearest = point
        if self.on_hover_callback is not None:
            self.on_hover_callback(nearest)
        return False

    def _on_leave(self, widget, event):
        if self.on_hover_callback is not None:
            self.on_hover_callback(None)
        return False
