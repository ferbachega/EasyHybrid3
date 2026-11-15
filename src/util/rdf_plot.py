#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: g(r) plot of the RDF Analysis window
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description (2026-10-21):
#      The RDF window used easyplot.XYPlot, made for energy profiles: its
#      axis margins put the g(r) axis below zero (-2 .. 5) and it has no
#      axis titles. This plot draws g(r) from 0 (left axis, solid) and the
#      running coordination number N(r) (right axis, dashed), with axis
#      titles, a hover readout and PNG export.
#

import math

import cairo
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk

_G_COLOR = (0.10, 0.10, 0.10)
_N_COLOR = (0.18, 0.42, 0.78)


def _nice_step(span, n=6):
    """ A 1/2/5 x 10^k step giving about n intervals over span. """
    if span <= 0: return 1.0
    raw = span / float(n)
    exp = math.floor(math.log10(raw))
    for m in (1, 2, 5, 10):
        if m * 10 ** exp >= raw:
            return m * 10 ** exp
    return 10 ** (exp + 1)


class RDFPlot(Gtk.DrawingArea):
    """ set_data(r, g, n=None, title=None); export_png(path). """

    def __init__(self):
        super().__init__()
        self.r, self.g, self.n, self.title = [], [], None, None
        self.hover = None
        self.add_events(Gdk.EventMask.POINTER_MOTION_MASK | Gdk.EventMask.LEAVE_NOTIFY_MASK)
        self.connect("draw", self.on_draw)
        self.connect("motion-notify-event", self._on_motion)
        self.connect("leave-notify-event", lambda w, e: self._set_hover(None))
        self.set_size_request(480, 300)

    def set_data(self, r, g, n=None, title=None):
        self.r, self.g = [float(x) for x in r], [float(y) for y in g]
        self.n = [float(y) for y in n] if n is not None else None
        self.title = title
        self.queue_draw()

    # . geometry
    def _frame(self, width, height):
        left, right = 64, (64 if self.n else 20)
        top, bottom = (34 if self.title else 16), 50
        return left, top, max(width - left - right, 10), max(height - top - bottom, 10)

    def _ranges(self):
        xmax = max(self.r) if self.r else 1.0
        gmax = max(self.g) if self.g else 1.0
        gstep = _nice_step(max(gmax, 1.0) * 1.08)
        gtop = math.ceil(max(gmax, 1.0) * 1.08 / gstep) * gstep
        ntop, nstep = None, None
        if self.n:
            nmax = max(self.n) or 1.0
            nstep = _nice_step(nmax * 1.05)
            ntop = math.ceil(nmax * 1.05 / nstep) * nstep
        return xmax, gtop, gstep, ntop, nstep

    def on_draw(self, widget, cr, width=None, height=None):
        width = width or widget.get_allocated_width()
        height = height or widget.get_allocated_height()
        cr.set_source_rgb(1, 1, 1); cr.paint()
        if not self.r:
            return False
        x0, y0, w, h = self._frame(width, height)
        xmax, gtop, gstep, ntop, nstep = self._ranges()
        X = lambda x: x0 + x / xmax * w
        Yg = lambda y: y0 + h - y / gtop * h
        cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        cr.set_font_size(11)
        # . grid and ticks
        cr.set_line_width(1)
        xstep = _nice_step(xmax, 8)
        k = 0
        while k * xstep <= xmax + 1e-9:
            x = X(k * xstep)
            cr.set_source_rgb(0.9, 0.9, 0.9); cr.move_to(x, y0); cr.line_to(x, y0 + h); cr.stroke()
            cr.set_source_rgb(0.2, 0.2, 0.2)
            text = "{:g}".format(round(k * xstep, 6)); e = cr.text_extents(text)
            cr.move_to(x - e.width / 2, y0 + h + 16); cr.show_text(text)
            k += 1
        k = 0
        while k * gstep <= gtop + 1e-9:
            y = Yg(k * gstep)
            cr.set_source_rgb(0.9, 0.9, 0.9); cr.move_to(x0, y); cr.line_to(x0 + w, y); cr.stroke()
            cr.set_source_rgb(0.2, 0.2, 0.2)
            text = "{:g}".format(round(k * gstep, 6)); e = cr.text_extents(text)
            cr.move_to(x0 - e.width - 6, y + e.height / 2); cr.show_text(text)
            k += 1
        # . g(r) = 1 reference
        cr.set_source_rgb(0.6, 0.6, 0.6); cr.set_dash([4, 4]); cr.move_to(x0, Yg(1)); cr.line_to(x0 + w, Yg(1)); cr.stroke(); cr.set_dash([])
        # . curves
        cr.save(); cr.rectangle(x0, y0, w, h); cr.clip()
        if self.n:
            Yn = lambda y: y0 + h - y / ntop * h
            cr.set_source_rgb(*_N_COLOR); cr.set_line_width(1.6); cr.set_dash([6, 4])
            for i, (x, y) in enumerate(zip(self.r, self.n)):
                (cr.move_to if i == 0 else cr.line_to)(X(x), Yn(y))
            cr.stroke(); cr.set_dash([])
        cr.set_source_rgb(*_G_COLOR); cr.set_line_width(2.0)
        for i, (x, y) in enumerate(zip(self.r, self.g)):
            (cr.move_to if i == 0 else cr.line_to)(X(x), Yg(y))
        cr.stroke(); cr.restore()
        # . frame
        cr.set_source_rgb(0.2, 0.2, 0.2); cr.set_line_width(1.2); cr.rectangle(x0, y0, w, h); cr.stroke()
        # . right axis, N(r)
        if self.n:
            k = 0
            while k * nstep <= ntop + 1e-9:
                y = y0 + h - k * nstep / ntop * h
                cr.set_source_rgb(*_N_COLOR)
                text = "{:g}".format(round(k * nstep, 6)); e = cr.text_extents(text)
                cr.move_to(x0 + w + 6, y + e.height / 2); cr.show_text(text)
                k += 1
        # . axis titles
        cr.set_font_size(12)
        cr.set_source_rgb(0.1, 0.1, 0.1)
        text = "r (Å)"; e = cr.text_extents(text); cr.move_to(x0 + w / 2 - e.width / 2, y0 + h + 38); cr.show_text(text)
        cr.save(); text = "g(r)"; e = cr.text_extents(text); cr.move_to(18, y0 + h / 2 + e.width / 2); cr.rotate(-math.pi / 2); cr.show_text(text); cr.restore()
        if self.n:
            cr.save(); cr.set_source_rgb(*_N_COLOR); text = "N(r)  (dashed)"; e = cr.text_extents(text)
            cr.move_to(x0 + w + 50, y0 + h / 2 + e.width / 2); cr.rotate(-math.pi / 2); cr.show_text(text); cr.restore()
        if self.title:
            cr.set_font_size(13); cr.set_source_rgb(0.1, 0.1, 0.1)
            e = cr.text_extents(self.title); cr.move_to(x0 + w / 2 - e.width / 2, y0 - 12); cr.show_text(self.title)
        # . hover readout
        if self.hover is not None:
            i = self.hover
            text = "r = {:.2f} Å   g = {:.3f}".format(self.r[i], self.g[i]) + ("   N = {:.2f}".format(self.n[i]) if self.n else "")
            cr.set_font_size(11); e = cr.text_extents(text)
            bx, by = x0 + w - e.width - 14, y0 + 8
            cr.set_source_rgba(1, 1, 1, 0.9); cr.rectangle(bx - 4, by - 2, e.width + 10, e.height + 8); cr.fill()
            cr.set_source_rgb(0.1, 0.1, 0.1); cr.move_to(bx, by + e.height + 1); cr.show_text(text)
            cr.set_source_rgb(0.85, 0.2, 0.2); cr.arc(X(self.r[i]), Yg(self.g[i]), 3.5, 0, 2 * math.pi); cr.fill()
        return False

    def _set_hover(self, index):
        if index != self.hover:
            self.hover = index; self.queue_draw()
        return False

    def _on_motion(self, widget, event):
        if not self.r: return False
        x0, y0, w, h = self._frame(widget.get_allocated_width(), widget.get_allocated_height())
        if not (x0 <= event.x <= x0 + w):
            return self._set_hover(None)
        x = (event.x - x0) / w * max(self.r)
        return self._set_hover(min(range(len(self.r)), key=lambda i: abs(self.r[i] - x)))

    def export_png(self, path, width=900, height=540):
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        cr = cairo.Context(surface)
        hover, self.hover = self.hover, None
        self.on_draw(self, cr, width, height)
        self.hover = hover
        surface.write_to_png(path)
