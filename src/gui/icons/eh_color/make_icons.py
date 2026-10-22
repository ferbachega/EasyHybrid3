#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid toolbar icons, "assorted colours" variant of icons/eh_style.
#
#  2026-10-22. Same drawings and outline style as eh_style (flat fills, thick
#  dark outline), but every icon has its own colours, close to the original
#  EasyHybrid icons (icons/icon_*.png): red/blue/cyan atoms, blue geometry
#  well, red dynamics well, purple umbrella, brown NEB, red E, blue Psi, ...
#
#  Drawn on a 24 x 24 grid and written as SVG at 32 x 32 px (svg/, the source
#  for future edits) and PNG at 32 px (png/32/, what MainWindow.glade loads:
#  PNG, not SVG, so the toolbar does not depend on the librsvg pixbuf loader).
#
#      python3 make_icons.py              # svg/ + png/32/
#      python3 make_icons.py --preview    # also preview_light.png / preview_dark.png
#
import math
import os
import sys

import cairo

HERE = os.path.dirname(os.path.abspath(__file__))
SVG_SIZE = 32
PNG_SIZES = (32,)

# . palette
INK     = (0.08, 0.08, 0.08)
PAPER   = (0.965, 0.961, 0.937)
WHITE   = (1.0, 1.0, 1.0)
GREY    = (0.62, 0.62, 0.62)
DGREY   = (0.30, 0.30, 0.30)
RED     = (0.890, 0.200, 0.180)
PINK    = (0.973, 0.643, 0.643)
ORANGE  = (0.976, 0.600, 0.180)
YELLOW  = (0.988, 0.835, 0.310)
LYELLOW = (0.988, 0.937, 0.620)
GREEN   = (0.290, 0.700, 0.270)
LGREEN  = (0.600, 0.850, 0.450)
CYAN    = (0.150, 0.720, 0.900)
BLUE    = (0.200, 0.420, 0.780)
LBLUE   = (0.690, 0.800, 0.930)
PURPLE  = (0.560, 0.250, 0.650)
LPURPLE = (0.800, 0.640, 0.860)
BROWN   = (0.560, 0.330, 0.150)
LBROWN  = (0.850, 0.640, 0.420)
SKIN    = (0.973, 0.847, 0.690)
LW      = 1.5


# ----------------------------------------------------------------------------
#  primitives
# ----------------------------------------------------------------------------
def fill_stroke(cr, fill, lw=LW):
    if fill is not None:
        cr.set_source_rgb(*fill); cr.fill_preserve()
    cr.set_source_rgb(*INK); cr.set_line_width(lw); cr.stroke()


def atom(cr, x, y, r=2.6, color=GREY, lw=1.2):
    cr.new_sub_path(); cr.arc(x, y, r, 0, 2 * math.pi)
    fill_stroke(cr, color, lw)


def bond(cr, x1, y1, x2, y2, lw=1.4, color=INK):
    cr.set_source_rgb(*color); cr.set_line_width(lw)
    cr.move_to(x1, y1); cr.line_to(x2, y2); cr.stroke()


def rrect(cr, x, y, w, h, r=2.0):
    cr.new_sub_path()
    cr.arc(x + w - r, y + r, r, -math.pi / 2, 0)
    cr.arc(x + w - r, y + h - r, r, 0, math.pi / 2)
    cr.arc(x + r, y + h - r, r, math.pi / 2, math.pi)
    cr.arc(x + r, y + r, r, math.pi, 3 * math.pi / 2)
    cr.close_path()


def poly(cr, pts, fill, close=True, lw=LW):
    cr.move_to(*pts[0])
    for p in pts[1:]: cr.line_to(*p)
    if close: cr.close_path()
    fill_stroke(cr, fill, lw)


def text_glyph(cr, text, x, y, size, fill, weight=cairo.FONT_WEIGHT_BOLD, family="DejaVu Sans", lw=1.0):
    cr.select_font_face(family, cairo.FONT_SLANT_NORMAL, weight)
    cr.set_font_size(size)
    e = cr.text_extents(text)
    cr.move_to(x - e.width / 2 - e.x_bearing, y - e.height / 2 - e.y_bearing)
    cr.text_path(text)
    fill_stroke(cr, fill, lw)


def molecule(cr, ox, oy, s=1.0, colors=(GREY, RED, GREY)):
    a, b, c = (ox, oy + 4 * s), (ox + 4 * s, oy), (ox + 8 * s, oy + 3 * s)
    bond(cr, *a, *b); bond(cr, *b, *c)
    for p, col in zip((a, b, c), colors): atom(cr, *p, 1.8 * s, col)


def document(cr, x=4, y=2, w=14, h=19, fold=4, fold_color=LBLUE):
    poly(cr, [(x, y), (x + w - fold, y), (x + w, y + fold), (x + w, y + h), (x, y + h)], WHITE)
    poly(cr, [(x + w - fold, y), (x + w - fold, y + fold), (x + w, y + fold)], fold_color, lw=1.1)


def arrow_down(cr, cx, top, bottom, half=4.0, shaft=1.6, fill=BLUE):
    hy = bottom - half - 0.5
    poly(cr, [(cx - shaft, top), (cx + shaft, top), (cx + shaft, hy), (cx + half, hy), (cx, bottom),
              (cx - half, hy), (cx - shaft, hy)], fill, lw=1.2)


def plus_badge(cr, cx, cy, r=4.6, fill=GREEN):
    a, t = r, r * 0.38
    poly(cr, [(cx - t, cy - a), (cx + t, cy - a), (cx + t, cy - t), (cx + a, cy - t), (cx + a, cy + t), (cx + t, cy + t),
              (cx + t, cy + a), (cx - t, cy + a), (cx - t, cy + t), (cx - a, cy + t), (cx - a, cy - t), (cx - t, cy - t)], fill, lw=1.2)


def plot_frame(cr, fill=WHITE):
    rrect(cr, 2, 2.5, 20, 19, 1.8); fill_stroke(cr, fill)


def curve(cr, pts, fill=None, close_to=None, lw=1.4):
    cr.move_to(*pts[0])
    for i in range(len(pts) - 1):
        p0 = pts[max(i - 1, 0)]; p1 = pts[i]; p2 = pts[i + 1]; p3 = pts[min(i + 2, len(pts) - 1)]
        c1 = (p1[0] + (p2[0] - p0[0]) / 6, p1[1] + (p2[1] - p0[1]) / 6)
        c2 = (p2[0] - (p3[0] - p1[0]) / 6, p2[1] - (p3[1] - p1[1]) / 6)
        cr.curve_to(*c1, *c2, *p2)
    if close_to is not None:
        cr.line_to(pts[-1][0], close_to); cr.line_to(pts[0][0], close_to); cr.close_path()
        fill_stroke(cr, fill, lw)
    else:
        cr.set_source_rgb(*INK); cr.set_line_width(lw); cr.stroke()


def pencil(cr, x0, y0, x1, y1, w=3.2, body=ORANGE):
    dx, dy = x1 - x0, y1 - y0; L = math.hypot(dx, dy); ux, uy = dx / L, dy / L; nx, ny = -uy, ux
    tip = 4.0
    bx, by = x0 + ux * tip, y0 + uy * tip
    poly(cr, [(bx + nx * w, by + ny * w), (x1 + nx * w, y1 + ny * w), (x1 - nx * w, y1 - ny * w), (bx - nx * w, by - ny * w)], body, lw=1.2)
    poly(cr, [(x0, y0), (bx + nx * w, by + ny * w), (bx - nx * w, by - ny * w)], SKIN, lw=1.2)
    ex, ey = x1 - ux * 2.6, y1 - uy * 2.6
    poly(cr, [(ex + nx * w, ey + ny * w), (x1 + nx * w, y1 + ny * w), (x1 - nx * w, y1 - ny * w), (ex - nx * w, ey - ny * w)], PINK, lw=1.2)


def curved_arrow(cr, color, mirror=False):
    cr.save()
    if mirror:
        cr.translate(24, 0); cr.scale(-1, 1)
    cr.set_line_cap(cairo.LINE_CAP_ROUND)
    for width, c in ((6.0, INK), (3.4, color)):
        cr.new_path(); cr.arc_negative(12.5, 13.5, 6.5, math.radians(40), math.radians(-178))
        cr.set_source_rgb(*c); cr.set_line_width(width); cr.stroke()
    cr.set_line_cap(cairo.LINE_CAP_BUTT)
    poly(cr, [(1.5, 12.0), (10.5, 12.0), (6.0, 19.5)], color, lw=1.3)
    cr.restore()


# ----------------------------------------------------------------------------
#  icons (24 x 24 grid)
# ----------------------------------------------------------------------------
def icon_new_project(cr):
    document(cr); molecule(cr, 7, 7.5)
    plus_badge(cr, 18, 18, fill=GREEN)

def icon_open(cr):
    poly(cr, [(2, 6), (8.5, 6), (10.5, 8), (21, 8), (21, 20), (2, 20)], YELLOW)
    poly(cr, [(2, 20), (5, 11), (23, 11), (20, 20)], LYELLOW, lw=1.3)
    molecule(cr, 9, 12.5, 0.8)

def icon_save(cr):
    document(cr, 3, 2, 13, 17); molecule(cr, 5.5, 6.5, 0.85)
    arrow_down(cr, 18, 9, 22, half=4.2, fill=BLUE)

def icon_save_as(cr):
    document(cr, 6, 1.5, 12, 15, 3.5)
    document(cr, 3, 4, 12, 16, 3.5); molecule(cr, 5.5, 8.5, 0.8)
    arrow_down(cr, 18.5, 10, 22.5, half=4.2, fill=ORANGE)

def icon_terminal(cr):
    rrect(cr, 2, 3.5, 20, 17, 2); fill_stroke(cr, (0.16, 0.17, 0.20))
    cr.set_source_rgb(*WHITE); cr.set_line_width(2.0); cr.set_line_cap(cairo.LINE_CAP_ROUND); cr.set_line_join(cairo.LINE_JOIN_ROUND)
    cr.move_to(6, 9); cr.line_to(9.5, 12); cr.line_to(6, 15); cr.stroke()
    cr.set_source_rgb(*CYAN)
    cr.move_to(11.5, 15.5); cr.line_to(17, 15.5); cr.stroke()
    cr.set_line_cap(cairo.LINE_CAP_BUTT)

def icon_play(cr):
    poly(cr, [(6, 3.5), (20, 12), (6, 20.5)], BLUE)

def icon_edit_build(cr):
    pts = [(9 + 6 * math.cos(math.radians(a)), 13 + 6 * math.sin(math.radians(a))) for a in range(-90, 270, 60)]
    for i in range(6): bond(cr, *pts[i], *pts[(i + 1) % 6])
    for k, p in enumerate(pts): atom(cr, *p, 1.9, RED if k == 0 else (BLUE if k == 3 else GREY))
    pencil(cr, 12.5, 15.5, 21.5, 2.5, w=2.4, body=YELLOW)

def icon_selections(cr):
    for col, y in zip((CYAN, ORANGE, CYAN), (5, 12, 19)):
        atom(cr, 5, y, 2.6, col)
        rrect(cr, 10, y - 1.6, 12, 3.2, 1.2); fill_stroke(cr, LBLUE if col is CYAN else LYELLOW, 1.1)

def icon_go_to_atom(cr):
    bond(cr, 4, 16, 7, 7); bond(cr, 7, 7, 15, 5); bond(cr, 4, 16, 10, 20)
    atom(cr, 4, 16, 2.2, GREY); atom(cr, 10, 20, 2.2, GREY); atom(cr, 15, 5, 2.2, GREY)
    atom(cr, 7, 7, 3.0, CYAN)
    poly(cr, [(10, 10), (21, 15), (16.5, 16.3), (19.5, 21), (17.2, 22.3), (14.3, 17.6), (11.2, 20.5)], RED, lw=1.3)

def icon_extend(cr):
    cr.set_dash([2.0, 1.6]); cr.new_sub_path(); cr.arc(12, 12, 10, 0, 2 * math.pi)
    cr.set_source_rgb(*BLUE); cr.set_line_width(1.4); cr.stroke(); cr.set_dash([])
    pts = [(6, 15), (10, 9), (15, 12), (18, 7)]
    for a, b in zip(pts, pts[1:]): bond(cr, *a, *b)
    for k, p in enumerate(pts): atom(cr, *p, 2.3, CYAN if k in (1, 2) else GREY)

def icon_qc_setup(cr):
    text_glyph(cr, "Ψ", 12, 12.5, 21, CYAN, family="DejaVu Serif", lw=1.0)

def icon_energy(cr):
    rrect(cr, 2.5, 2.5, 19, 19, 3); fill_stroke(cr, WHITE)
    text_glyph(cr, "E", 12, 12.3, 15, RED, family="DejaVu Serif", lw=1.0)

def _well(cr, fill, node, pts, node_xy):
    plot_frame(cr)
    curve(cr, pts, fill, close_to=21.5)
    if node: atom(cr, *node_xy, 2.4, node)

def icon_geometry_optimization(cr):
    _well(cr, BLUE, YELLOW, [(2, 5), (7, 13), (12, 17), (17, 13), (22, 5)], (12, 14.4))

def icon_molecular_dynamics(cr):
    plot_frame(cr)
    curve(cr, [(2, 9), (6.5, 16), (11, 12.5), (15.5, 16.5), (22, 8)], RED, close_to=21.5)
    cr.set_source_rgb(*DGREY); cr.set_line_width(1.2); cr.set_line_cap(cairo.LINE_CAP_ROUND)
    for y in (5.2, 7.4, 9.6):
        cr.move_to(4.5, y); cr.line_to(8.0, y); cr.stroke()
    cr.set_line_cap(cairo.LINE_CAP_BUTT)
    atom(cr, 12, 7.4, 2.8, YELLOW)

def icon_normal_modes(cr):
    zig = [(7.5, 12)] + [(8.8 + i * 1.5, 12 + (2.2 if i % 2 == 0 else -2.2)) for i in range(5)] + [(16.5, 12)]
    cr.set_source_rgb(*INK); cr.set_line_width(1.2); cr.move_to(*zig[0])
    for p in zig[1:]: cr.line_to(*p)
    cr.stroke()
    atom(cr, 5, 12, 3.2, ORANGE); atom(cr, 19, 12, 3.2, ORANGE)
    for x0 in (5, 19):
        cr.set_source_rgb(*RED); cr.set_line_width(1.3)
        cr.move_to(x0 - 3, 5); cr.line_to(x0 + 3, 5); cr.stroke()
        for xa, d in ((x0 - 3, 1), (x0 + 3, -1)):
            cr.move_to(xa + d * 1.6, 3.6); cr.line_to(xa, 5); cr.line_to(xa + d * 1.6, 6.4); cr.stroke()

def icon_scan(cr):
    plot_frame(cr)
    pts = [(2, 17), (7, 15), (12, 6.5), (17, 14), (22, 16)]
    curve(cr, pts, GREEN, close_to=21.5)
    for x, y in ((5.5, 15.6), (12, 6.6), (18.5, 14.6)): atom(cr, x, y - 0.4, 1.7, RED)

def icon_umbrella(cr):
    plot_frame(cr)
    for c, col in ((6.5, LPURPLE), (12, PURPLE), (17.5, LPURPLE)):
        curve(cr, [(c - 5, 21.5), (c - 2.6, 17), (c, 9), (c + 2.6, 17), (c + 5, 21.5)], col, close_to=21.5, lw=1.1)

def icon_neb(cr):
    plot_frame(cr)
    pts = [(2, 11), (5.5, 17), (10, 12), (12, 8), (14, 12), (18.5, 17), (22, 11)]
    curve(cr, pts, LBROWN, close_to=21.5)
    for x, y in ((5.5, 16.5), (9, 13), (12, 8.3), (15, 13), (18.5, 16.5)): atom(cr, x, y - 0.5, 1.6, BROWN)

def icon_export_image(cr):
    rrect(cr, 2, 3, 17, 14, 1.5); fill_stroke(cr, LBLUE)
    poly(cr, [(3, 16), (8, 9.5), (11.5, 13), (13.5, 11), (18, 16)], GREEN, lw=1.1)
    atom(cr, 14.5, 7, 1.7, YELLOW, lw=1.0)
    arrow_down(cr, 18.5, 11.5, 23, half=4.0, fill=BLUE)

def icon_process_manager(cr):
    rrect(cr, 2, 3, 20, 18, 2); fill_stroke(cr, WHITE)
    for y, frac, col in ((8, 1.0, GREEN), (12.5, 0.65, CYAN), (17, 0.3, ORANGE)):
        rrect(cr, 5, y - 1.5, 14, 3, 1.2); fill_stroke(cr, (0.90, 0.90, 0.88), 0.9)
        rrect(cr, 5, y - 1.5, 14 * frac, 3, 1.2); fill_stroke(cr, col, 0.9)

def icon_editing(cr):
    pencil(cr, 4, 20, 19.5, 4.5, w=3.0, body=ORANGE)

def icon_add(cr):
    plus_badge(cr, 12, 12, r=9.0, fill=GREEN)

def icon_delete(cr):
    t, a = 2.6, 8.5
    cr.save(); cr.translate(12, 12); cr.rotate(math.radians(45))
    poly(cr, [(-t, -a), (t, -a), (t, -t), (a, -t), (a, t), (t, t), (t, a), (-t, a), (-t, t), (-a, t), (-a, -t), (-t, -t)], RED)
    cr.restore()

def icon_transform(cr):
    for dx, dy in ((0, -1), (0, 1), (-1, 0), (1, 0)):
        tip = (12 + dx * 10, 12 + dy * 10); base = (12 + dx * 6.5, 12 + dy * 6.5)
        nx, ny = -dy, dx
        poly(cr, [tip, (base[0] + nx * 3, base[1] + ny * 3), (base[0] - nx * 3, base[1] - ny * 3)], LBLUE, lw=1.1)
    bond(cr, 9, 14, 12, 9.5); bond(cr, 12, 9.5, 15.5, 13.5)
    atom(cr, 9, 14, 2.3, GREY); atom(cr, 12, 9.5, 2.3, RED); atom(cr, 15.5, 13.5, 2.3, GREY)

def icon_undo(cr):
    curved_arrow(cr, LYELLOW)

def icon_redo(cr):
    curved_arrow(cr, PINK, mirror=True)

def icon_add_hydrogens(cr):
    bond(cr, 6, 15, 12, 11)
    atom(cr, 6, 15, 3.4, DGREY)
    atom(cr, 12.5, 10.5, 4.6, WHITE)
    text_glyph(cr, "H", 12.5, 10.6, 6.8, INK, lw=0.1)
    plus_badge(cr, 18.5, 18.5, r=4.3, fill=GREEN)

def icon_force_field_settings(cr):
    for y, x, col in ((5, 15, RED), (12, 7.5, BLUE), (19, 13, ORANGE)):
        rrect(cr, 2.5, y - 1.1, 19, 2.2, 1.0); fill_stroke(cr, (0.90, 0.90, 0.88), 1.0)
        atom(cr, x, y, 2.8, col)


def icon_bond_order(cr):
    # . two atoms joined by a double bond: one solid line, one dashed
    cr.save(); cr.translate(12, 12); cr.rotate(math.radians(-45))
    cr.set_source_rgb(*INK); cr.set_line_width(1.6)
    cr.move_to(-8, -1.8); cr.line_to(8, -1.8); cr.stroke()
    cr.set_dash([1.8, 1.4]); cr.move_to(-8, 1.8); cr.line_to(8, 1.8); cr.stroke(); cr.set_dash([])
    cr.restore()
    atom(cr, 5.5, 18.5, 4.2, BLUE); atom(cr, 18.5, 5.5, 4.2, RED)

def icon_dihedral(cr):
    a, b, c, d = (5, 4.5), (8, 12.5), (17, 13.5), (19.5, 21)
    poly(cr, [a, b, c], LBLUE, lw=1.0)
    poly(cr, [b, c, (20, 6)], LYELLOW, lw=1.0)
    for p, q in ((a, b), (b, c), (c, d)): bond(cr, *p, *q, lw=1.6)
    cr.set_source_rgb(*INK); cr.set_line_width(1.0); cr.set_dash([1.4, 1.2])
    cr.move_to(*c); cr.line_to(20, 6); cr.stroke(); cr.set_dash([])
    atom(cr, *a, 2.3, GREY); atom(cr, *b, 2.6, ORANGE); atom(cr, *c, 2.6, ORANGE); atom(cr, *d, 2.3, GREY)
    atom(cr, 20, 6, 1.6, WHITE, lw=1.0)

def icon_move(cr):
    atom(cr, 9, 9, 6.5, CYAN)
    # . pointing hand: index finger up-left onto the atom, folded fingers below
    poly(cr, [(10.5, 10.5), (12.6, 9.6), (15.8, 15.0), (17.2, 13.9), (18.6, 14.7), (19.8, 14.0), (21.4, 15.2),
              (22.4, 15.0), (23.2, 18.5), (21.6, 22.5), (16.0, 23.0), (12.8, 18.4), (11.2, 16.0), (12.4, 15.0), (14.0, 16.0)],
         SKIN, lw=1.2)
    cr.set_source_rgb(*INK); cr.set_line_width(0.9)
    cr.move_to(17.2, 13.9); cr.line_to(17.9, 16.6); cr.stroke()
    cr.move_to(19.8, 14.0); cr.line_to(20.3, 16.6); cr.stroke()

def icon_optimize(cr):
    _well(cr, GREEN, YELLOW, [(2, 4), (7, 13), (12, 17), (17, 13), (22, 4)], (12, 14.4))

def icon_solvate(cr):
    # . a molecule inside a box of water
    poly(cr, [(4, 8), (9, 3), (21, 3), (21, 15), (16, 20), (4, 20)], LBLUE)
    poly(cr, [(4, 8), (9, 3), (21, 3), (16, 8)], (0.82, 0.89, 0.97), lw=1.1)
    poly(cr, [(16, 8), (21, 3), (21, 15), (16, 20)], (0.55, 0.70, 0.88), lw=1.1)
    rrect(cr, 4, 8, 12, 12, 0.01); fill_stroke(cr, None, 1.1)
    for x, y in ((6.5, 11), (13.5, 17.5), (6.5, 17.5), (13.5, 11)):
        atom(cr, x, y, 1.1, RED, lw=0.8)
    bond(cr, 8.3, 15.5, 11.7, 12.5, lw=1.2)
    atom(cr, 8.3, 15.5, 2.0, GREY); atom(cr, 11.7, 12.5, 2.0, ORANGE)


ICONS = [
    # . main toolbar
    ("new_project", icon_new_project), ("open", icon_open), ("save", icon_save), ("save_as", icon_save_as),
    ("terminal", icon_terminal), ("trajectory_player", icon_play), ("edit_and_build", icon_edit_build),
    ("selections", icon_selections), ("go_to_atom", icon_go_to_atom), ("extend_selection", icon_extend),
    ("qc_setup", icon_qc_setup), ("energy", icon_energy), ("geometry_optimization", icon_geometry_optimization),
    ("molecular_dynamics", icon_molecular_dynamics), ("normal_modes", icon_normal_modes), ("pes_scan", icon_scan),
    ("umbrella_sampling", icon_umbrella), ("neb", icon_neb), ("export_image", icon_export_image),
    ("process_manager", icon_process_manager),
    # . Builder toolbar
    ("editing", icon_editing), ("add", icon_add), ("delete", icon_delete), ("transform", icon_transform),
    ("undo", icon_undo), ("redo", icon_redo), ("add_hydrogens", icon_add_hydrogens),
    ("force_field_settings", icon_force_field_settings), ("bond_order", icon_bond_order),
    ("move", icon_move), ("dihedral", icon_dihedral), ("optimize", icon_optimize), ("solvate", icon_solvate),
]


def render(func, size):
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, size, size)
    cr = cairo.Context(surface); cr.scale(size / 24.0, size / 24.0)
    cr.set_line_join(cairo.LINE_JOIN_ROUND)
    func(cr)
    return surface


def preview(path, background, size=32, per_row=14):
    pad = 10; rows = (len(ICONS) + per_row - 1) // per_row
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, pad + per_row * (size + pad), pad + rows * (size + pad))
    cr = cairo.Context(surface); cr.set_source_rgb(*background); cr.paint()
    for k, (name, func) in enumerate(ICONS):
        cr.set_source_surface(render(func, size), pad + (k % per_row) * (size + pad), pad + (k // per_row) * (size + pad))
        cr.paint()
    surface.write_to_png(path)


def main():
    os.makedirs(os.path.join(HERE, "svg"), exist_ok=True)
    for name, func in ICONS:
        svg = cairo.SVGSurface(os.path.join(HERE, "svg", name + ".svg"), SVG_SIZE, SVG_SIZE)
        svg.set_document_unit(cairo.SVGUnit.PX)
        cr = cairo.Context(svg); cr.scale(SVG_SIZE / 24.0, SVG_SIZE / 24.0)
        cr.set_line_join(cairo.LINE_JOIN_ROUND); func(cr); svg.finish()
        for size in PNG_SIZES:
            folder = os.path.join(HERE, "png", str(size)); os.makedirs(folder, exist_ok=True)
            render(func, size).write_to_png(os.path.join(folder, "icon_{}.png".format(name)))
    if "--preview" in sys.argv:
        preview(os.path.join(HERE, "preview_light.png"), (0.95, 0.95, 0.95))
        preview(os.path.join(HERE, "preview_dark.png"), (0.20, 0.20, 0.20))
    print("{} icons written to {}".format(len(ICONS), HERE))


if __name__ == "__main__":
    main()
