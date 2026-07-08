#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Copyright 2022-2025 Fernando Bachega
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
#      Provides functions for selecting atoms and residues in pDynamo systems
#      to facilitate QM/MM partitioning and molecular simulations.
#
#===============================================================================
#  SESSION CHANGELOG -- Surface rendering overhaul (chat-assisted session)
#===============================================================================
#  Everything below documents a single, continuous round of changes made to
#  this file (and, in a few cases, to representations.py / vismol_glcore.py /
#  shaders/surface.py in the graphics_engine submodule, and to a new file,
#  src/util/cube_reader.py) during one collaborative debugging/feature session.
#  Comments are kept in English here regardless of the surrounding code's
#  language, at the author's request, to make the history easier to follow
#  for anyone reading the module later. Existing Portuguese comments sprinkled
#  through the functions below were written inline as each fix/feature landed;
#  this block is the consolidated, chronological summary.
#
#  1. TRIANGLE-FILLED SURFACE RENDERING (was wireframe-only)
#     ----------------------------------------------------------------------
#     Originally, SurfaceRepresentation (in representations.py) rendered every
#     surface (orbitals, density, potential) as a GL_LINES wireframe, even
#     though a complete Phong-lit triangle shader already existed in
#     shaders/surface.py but was never wired up (the "surface" shader_program
#     key was bound to the *lines* shader pair in vismol_glcore.py, and the
#     draw call used GL_LINES instead of GL_TRIANGLES).
#     Fixed by: switching the shader_program binding to the triangle shaders,
#     switching the draw call to GL_TRIANGLES, and -- because the geometry
#     shader depended on a "vert_normal" vertex attribute that had no backing
#     VBO at the time -- making the geometry shader compute a flat per-face
#     normal from the triangle's own edges (cross product) as an interim fix.
#     (This file was not touched for this specific step; see representations.py
#     and shaders/surface.py in the graphics_engine submodule.)
#
#  2. RENDER MODE TOGGLE: filled surface vs. wireframe (user-selectable)
#     ----------------------------------------------------------------------
#     Added a "Wireframe" checkbutton (self.chk_surface_wireframe) next to the
#     surface-type combobox. Rather than maintaining two separate shader/VAO
#     pipelines, both modes reuse the exact same triangulated mesh and shader;
#     what changes is only the OpenGL rasterization mode, toggled via
#     glPolygonMode(GL_FRONT_AND_BACK, GL_FILL | GL_LINE) inside
#     SurfaceRepresentation.draw_representation() (representations.py), driven
#     by a new self.render_mode attribute and set_render_mode() method on that
#     class. Handler here: on_surface_wireframe_toggled(), which walks every
#     vismol object flagged is_surface == True in self.vm_session.vm_objects_dic
#     and calls set_render_mode() on each of its representations.
#
#  3. TRANSPARENCY / OPACITY SLIDER
#     ----------------------------------------------------------------------
#     Added an "Opacity" Gtk.Scale (0-100%, default 100% = fully opaque, i.e.
#     unchanged behaviour) next to the wireframe checkbutton. GL_BLEND and the
#     corresponding blend function were already enabled in
#     draw_representation() from the start, but had no visible effect because
#     the fragment shader always wrote alpha = 1.0. Added a "surf_alpha"
#     uniform to fragment_shader_surface (shaders/surface.py), a matching
#     self.alpha attribute + set_alpha() method on SurfaceRepresentation, and
#     -- when alpha < 1.0 -- a temporary glDepthMask(GL_FALSE) around the draw
#     call (restored to GL_TRUE right after) to avoid the surface's own back
#     faces occluding its front faces in an ugly way. Handler here:
#     on_surface_opacity_changed(), same object-walking pattern as the
#     wireframe toggle.
#
#  4. SMOOTH SHADING (per-vertex normals) vs. FLAT SHADING (per-face normals)
#     ----------------------------------------------------------------------
#     Item 1 above left the renderer with flat shading only (visible facets
#     from the marching-cubes triangulation). Added a "Smooth shading"
#     checkbutton (self.chk_surface_smooth) to toggle genuine per-vertex
#     normal interpolation. This required:
#       - A per-vertex normal buffer that never previously existed (the
#         "vert_normal" attribute was declared in the vertex shader from the
#         start but had no data behind it -- see representations.py:
#         _make_gl_normal_buffer / _load_normal_vbo).
#       - geometry_shader_surface gained a "smooth_shading" uniform (int):
#         0 (default) keeps using the flat, per-face normal computed from the
#         triangle's own edges; non-zero switches to the interpolated
#         per-vertex normal instead.
#       - vertex_shader_surface now transforms vert_normal into view space
#         (mat3(view_mat * model_mat) * vert_normal) instead of passing it
#         through untransformed, so it stays consistent with the flat normal
#         (which is derived directly from view-space triangle edges).
#       - HERE, in this file: the actual per-vertex normal values are computed
#         inside surface_parser() / surface_parser_mep() below. First attempt
#         was a hand-written compute_smooth_normals() (pure-Python loop:
#         accumulate each triangle's face normal into its three vertices,
#         then normalize -- classic area-weighted vertex-normal averaging).
#         That function was LATER REMOVED once it was discovered that
#         pDynamo3 already ships the exact same algorithm, compiled, as
#         PolygonalSurface.MakeVertexNormalsFromPolygonalNormals() (source:
#         github.com/pdynamo/pDynamo3, pScientific/Surfaces/PolygonalSurface.py).
#         It is allocated but never populated by
#         QCGridPropertyGenerator.Isosurface() (which only ever calls
#         MakePolygonNormals()), so calling it ourselves, once, right after
#         getting the surface object, was a safe, correct, and much cheaper
#         substitute for the hand-rolled Python version. Both parser
#         functions now do exactly that.
#     Handler here: on_surface_smooth_toggled(), same object-walking pattern.
#
#  5. BUG FIXES found while implementing/exercising the "Electrostatic
#     Potential" surface type (index 1 in the combobox)
#     ----------------------------------------------------------------------
#     a) A stray `return False` sat as literally the first statement inside
#        the `elif index == 1:` branch of on_render_button() -- clicking
#        "Render" with "Electrostatic Potential" selected did *nothing*,
#        silently, because the function bailed out before building the job
#        list or calling generate_grid_parallel() at all. Removed.
#     b) In generate_grid_parallel()'s `elif _type == 'potential':` branch,
#        the positive-lobe isosurface was being coloured with `color_minus`
#        instead of `color_plus` (simple copy-paste slip from the orbital
#        branch above it). Fixed.
#     c) Same branch, the *negative*-lobe isosurface generation called
#        `generator.Isosurface(_OrbitalTag, ...)` -- `_OrbitalTag` being a
#        leftover variable name copied from the orbital branch, not the
#        string 'potential' that this branch actually needs. Fixed to use
#        'potential' explicitly.
#     d) A real pDynamo3 tag-collision bug: '_IsosurfaceTag' was set equal to
#        'potential' in the parameters dict built by on_render_button(). The
#        first call to generator.Isosurface('potential', isovalue,
#        tag=_IsosurfaceTag) reads the raw potential grid from the
#        'potential' key and *overwrites that same key* with the resulting
#        isosurface object (since tag == 'potential' too). The *second*
#        Isosurface() call (for the negative lobe) then tries to read the
#        raw grid from 'potential' again, but finds the isosurface object
#        left there by the first call instead -- raising pDynamo's own
#        "Invalid QC property for isosurface generation" error. Confirmed
#        live via the user's traceback. Fixed by giving the isosurface
#        output its own, separate tag ('Isosurface', mirroring the pattern
#        already used -- correctly -- by the orbital branch, where
#        _OrbitalTag and _IsosurfaceTag are always two different strings)
#        and hard-coding the raw-grid tag to 'potential' inside
#        generate_grid_parallel() instead of reusing the (now different)
#        _IsosurfaceTag variable for both purposes.
#     e) Same latent tag-collision bug also exists in the 'density' branch
#        (_IsosurfaceTag == 'density' there too) but is harmless today
#        because that branch only calls Isosurface() once. Left as a
#        documented risk, not fixed, since fixing it wasn't asked for and
#        touching it isn't currently load-bearing.
#     f) The 'density' branch was also colouring its (only) isosurface with
#        `color_minus` instead of `color_plus`. Fixed in the same pass as (b).
#
#  6. MEP (Molecular Electrostatic Potential mapped onto a density isosurface)
#     ----------------------------------------------------------------------
#     Distinct from the existing "Electrostatic Potential" surface type
#     (which draws an isosurface of the potential FIELD itself, at a fixed
#     value): true MEP takes the DENSITY isosurface's geometry and colours
#     each vertex by the potential VALUE interpolated at that point -- a
#     continuous red/white/blue-style colour gradient, the standard
#     representation used by Gaussian/GaussView, Avogadro, VMD, etc.
#     Added as a 5th combobox entry, "MEP (density + potential)". Building
#     blocks, all new in this file:
#       - build_potential_interpolator(potentialProperty): wraps a
#         scipy.spatial.cKDTree nearest-neighbour lookup over the pDynamo
#         potential grid's raw (gridPoints, gridValues) pair. Nearest-
#         neighbour, not trilinear, was used here because the pDynamo grid
#         object doesn't expose its own [i,j,k] indexing order to Python --
#         only flat parallel point/value arrays.
#       - mep_colormap(values, vmin, vmax, cmap_name, reverse, percentile):
#         maps a 1-D array of potential values to RGB. Went through two
#         implementations: first a hand-written linear red-white-blue
#         interpolation, then a version using matplotlib.colormaps with
#         mcolors.TwoSlopeNorm (better perceptual centring at zero), and
#         FINALLY (current version) switched to the project's own
#         COLOR_MAPS dictionary (src/util/colormaps.py, confirmed identical
#         to a file the user uploaded) via a small custom vectorised
#         interpolator, _colormap_lookup(), instead of depending on
#         matplotlib at all for this feature. Also gained: (i) a
#         "percentile" parameter (default 2.0) so that automatic vmin/vmax
#         are taken from the 2nd/98th percentile of the data rather than
#         the raw min/max -- this matters a lot in practice, see item 8
#         below; (ii) a "reverse" flag, because COLOR_MAPS's diverging maps
#         (coolwarm, vik, berlin...) run low-to-high as blue-to-red, the
#         opposite of the chemistry convention (red = negative/electron-
#         rich, blue = positive/electron-poor) this feature wants by
#         default.
#       - surface_parser_mep(surface, vertex_colors): a variant of
#         surface_parser() that takes a pre-computed per-original-vertex
#         RGB array instead of one flat iso_color repeated for the whole
#         mesh.
#     GUI additions for MEP: self.entry_mep_vmin / self.entry_mep_vmax
#     (optional manual override of the colour-scale bounds; empty = keep
#     the automatic percentile-based behaviour), and self.cbx_mep_cmap
#     (a combobox populated directly from COLOR_MAPS.keys(), so any
#     colormap added later to that dict shows up here automatically with
#     no further code changes). All three are shown only when "MEP" is the
#     active combobox entry (see surface_combobox_change()).
#
#  7. A REAL RENDERING BUG found *because of* MEP: colours were not being
#     read per vertex at all
#     ----------------------------------------------------------------------
#     After wiring up MEP's continuous colour gradient, the whole surface
#     rendered as a single, nearly uniform, washed-out grey -- not a
#     gradient. Root cause (in representations.py, not this file): the
#     colour VBO was created with
#     _make_gl_color_buffer(colors, shader_program, instances=True), which
#     sets glVertexAttribDivisor(attr, 1) -- meaning the "vert_color"
#     attribute advances once per OpenGL *instance*, not once per vertex.
#     But draw_representation() issues a plain glDrawElements() call (not
#     glDrawElementsInstanced()), so there is effectively only ever "one
#     instance" -- with the divisor active, the GPU always sampled
#     colors[0] for literally every vertex in the mesh, silently, for as
#     long as every vertex happened to share the same colour anyway
#     (true for plain orbital/potential/density surfaces, which use one
#     flat iso_color). MEP was the first case with genuinely different
#     colours per vertex, which is what finally made the bug visible.
#     Fixed by dropping instances=True from that call.
#
#  8. A REAL DATA/SCALING BUG found *because of* MEP: a single outlier
#     value silently wrecking the whole colour scale
#     ----------------------------------------------------------------------
#     Even after fixing item 7, the surface still looked almost uniformly
#     grey. Live debug prints (temporarily added, then removed) showed
#     potential values ranging from about -0.25 to +54010 across ~226
#     vertices -- one single vertex (very likely sitting extremely close
#     to a nucleus, where the electrostatic potential genuinely diverges
#     like 1/r) had a magnitude roughly five orders of magnitude larger
#     than every other vertex. Because mep_colormap() originally used the
#     raw min()/max() of the data to set its colour-scale limits, that one
#     outlier alone defined the whole scale -- every "normal" vertex ended
#     up within a fraction of a percent of the exact centre of the colour
#     range, i.e. visually indistinguishable near-uniform grey. Fixed by
#     switching the *default* (no explicit vmin/vmax given) behaviour to
#     use the 2nd/98th percentile instead of the true min/max (see item 6);
#     the rare extreme vertices simply get clamped/saturated to the most
#     extreme colour instead of dragging the whole scale toward the middle.
#
#  9. EXTERNAL CUBE FILE IMPORT (e.g. orbitals/density/ESP exported from
#     ORCA via its orca_plot utility)
#     ----------------------------------------------------------------------
#     The combobox already had a 4th entry, "External", and a
#     GtkFileChooserButton (btn_external_file) already existed in the
#     .glade file and was already shown/hidden for that entry -- but it had
#     never been wired to any handler at all; selecting "External" and
#     clicking Render did nothing. This session implemented it properly:
#       - NEW FILE: src/util/cube_reader.py -- a small, dependency-free
#         (only numpy) Gaussian Cube format reader. Class CubeGrid (title,
#         comment, atoms, origin, voxel_vectors, dims, values) plus
#         read_cube_file(path) plus a CubeFileError exception with
#         specific, actionable messages (negative NATOMS = unsupported
#         "multi-cube"/multi-orbital-per-file variant; value count
#         mismatch = truncated/corrupt file; etc). Verified against two
#         hand-generated synthetic .cube files (a two-Gaussian "density"
#         and a two-soft-Coulomb-centre "potential") and against a
#         deliberately truncated file, all producing the expected results.
#       - cube_to_pdynamo_surface(cube_grid, isovalue), HERE in this file:
#         builds a pDynamo3 RegularGrid (pScientific.Geometry3) via
#         RegularGrid.FromDimensionData([{bins, binSize, lower}, ...]) and
#         a RealArrayND (pScientific.Arrays, via Array.FromIterable() +
#         Reshape()) directly from the cube data, then calls
#         MarchingCubes_Isosurface3D(grid, dataND, isovalue) -- the exact
#         same compiled marching-cubes routine that
#         QCGridPropertyGenerator.Isosurface() already uses internally for
#         every other surface type (confirmed by reading pDynamo3's own
#         source, pSimulation/QCGridProperties.py and
#         pScientific/Surfaces/__extensions__/pyrex/
#         pScientific.Surfaces.MarchingCubes.pyx, cloned fresh from
#         github.com/pdynamo/pDynamo3 for this purpose). The returned
#         PolygonalSurface is therefore a drop-in match for
#         surface_parser()/surface_parser_mep() -- no adapter code needed
#         downstream at all. Only limitation: requires an axis-aligned
#         (diagonal) voxel matrix; raises ValueError otherwise (cube files
#         with sheared/rotated axes are rare but do exist for some
#         programs; not handled here).
#       - build_potential_interpolator_from_cube(cube_grid): the external-
#         file analogue of build_potential_interpolator(), but -- unlike
#         the opaque pDynamo grid object, which forced a nearest-neighbour
#         KDTree approach -- a parsed CubeGrid gives an explicit regular
#         grid structure (known origin + per-axis spacing), so this uses
#         scipy.interpolate.RegularGridInterpolator for genuine trilinear
#         interpolation instead, which is strictly more accurate.
#       - _generate_external_cube_surface(parameters) plus an early-return
#         guard at the very top of generate_grid_parallel(): every other
#         surface type unconditionally calls apply_coords_to_system(),
#         system.Energy(), and QCGridPropertyGenerator.FromSystem(system)
#         before even looking at which `_type` was requested -- all of
#         which require a live pDynamo QC system that simply does not
#         exist for an externally-supplied cube file (system/coords are
#         None in that job). The external-cube case now returns early,
#         before any of that QC setup runs, reading only from disk +
#         calling the marching-cubes adapter above. If a second, optional
#         potential cube file is also supplied, the density surface is
#         coloured via the same MEP machinery (item 6) instead of a flat
#         iso_color.
#       - GUI: reused the existing btn_external_file (a
#         GtkFileChooserButton from the .glade file -- NOT a plain
#         Gtk.Button; it emits "file-set", not "clicked" -- an early
#         attempt to wire it up as if it were a normal button crashed with
#         "TypeError: unknown signal name: clicked" the first time this
#         was actually run, live, in the user's environment). Added a
#         second, optional GtkFileChooserButton created purely in Python,
#         self.btn_external_potential_file, styled identically, for the
#         optional potential cube (MEP colouring for externally-imported
#         data too). Both widgets, plus their matching labels, are shown
#         only while "External" is the active combobox entry.
#
#  10. TWO GENUINE INITIALIZATION-ORDER BUGS, both only surfaced once the
#      window was actually opened in the user's live environment
#      ----------------------------------------------------------------------
#      a) self.cbx_surface_type.set_active(0), called right after the
#         combobox was populated, fires GTK's "changed" signal
#         SYNCHRONOUSLY and IMMEDIATELY -- i.e. surface_combobox_change()
#         ran before several widgets it references (label_mep_vmin,
#         cbx_mep_cmap, label_external_potential, ...) had even been
#         created yet further down in the same __init__ block, raising
#         AttributeError. Fixed by moving that set_active(0) call to the
#         very end of the widget-construction sequence in open_window(),
#         after every widget any combobox-change handler could possibly
#         touch already exists.
#      b) (documented under item 9 above) the "unknown signal name:
#         clicked" crash from treating a GtkFileChooserButton as if it
#         were a plain Gtk.Button.
#
#  11. MARCHING-CUBES "GHOST VERTEX" ARTIFACT (spurious long triangles
#      shooting off from the surface to a single distant point -- see the
#      screenshot the user attached)
#      ----------------------------------------------------------------------
#      First hypothesis (WRONG, or at least incomplete): a stray vertex
#      landing at exactly the *world* origin (0,0,0), caused by some
#      ambiguous marching-cubes cube configuration (inside
#      pDynamo3's MarchingCubes.c, function AddInteriorVertex, used for
#      cube cases 7/10/12/13 of the classic algorithm) leaving a vertex
#      slot untouched at its zero-initialised default. A first filter,
#      _is_degenerate_vertex(), discarded any triangle referencing a
#      vertex whose three coordinates were EXACTLY 0.0 (deliberately exact
#      equality, not a small-distance tolerance -- a real, continuously-
#      interpolated vertex essentially never lands on the exact bit
#      pattern of 0.0 in all three components simultaneously, even for a
#      molecule that happens to sit right at the coordinate origin, as our
#      synthetic test .cube does on purpose).
#      That first filter turned out to be insufficient in the user's real
#      case (still produced the artifact, per the screenshot). Re-reading
#      MarchingCubes.c more carefully (lines ~280-293) revealed the real
#      mechanism: raw vertex positions are computed in GRID-INDEX units
#      and only scaled (* binSize) and TRANSLATED (+ midPointLower -- the
#      grid's OWN origin corner, not the world origin) at the very end.
#      A never-written vertex slot (raw value (0,0,0)) therefore ends up
#      at the grid's own midPointLower corner after that final transform
#      -- which is essentially always somewhere far from the molecule,
#      but essentially NEVER at the world's (0,0,0) unless the grid
#      happens to have been defined with its own origin exactly there
#      (true, by construction, for the synthetic test cube -- which is
#      exactly why the first, position-based filter looked like it worked
#      during earlier testing, but doesn't generalise).
#      FINAL fix: switched from a position-based test to a purely
#      GEOMETRIC one. New _compute_valid_polygon_mask(polygons, vertices,
#      size_factor=8.0) computes the longest edge of every triangle in the
#      mesh, takes the MEDIAN of all of those lengths as a self-calibrating
#      estimate of "one normal marching-cubes triangle's size" (no need to
#      know the grid spacing explicitly), and discards any triangle whose
#      longest edge exceeds size_factor times that median -- a real
#      isosurface triangle is never much larger than one grid cell; a
#      triangle connecting a good vertex to a distant ghost vertex is,
#      by construction, enormous by comparison. Used in both
#      surface_parser() and surface_parser_mep(); _is_degenerate_vertex()
#      is kept in the file only as a documented, no-longer-called
#      reference for how the (incomplete) first attempt worked.
#      Verified against a hand-built mock mesh (27 normal small triangles
#      plus 3 artificially attached to one distant "ghost" vertex): the
#      mask correctly flags exactly those 3 as invalid and keeps the other
#      27.
#
#  12. SCIPY DEPENDENCY REMOVED ENTIRELY (both remaining uses replaced with
#      hand-written, pure-numpy equivalents)
#      ----------------------------------------------------------------------
#      Triggered by two things happening back to back: (a) a second machine
#      hit `ModuleNotFoundError: No module named 'scipy'` -- and because
#      this file is imported unconditionally very early during EasyHybrid
#      startup, that took the WHOLE application down, not just the surface
#      window (a "soft"/optional-import version of scipy was tried first,
#      see the now-superseded SCIPY_AVAILABLE / _require_scipy() approach
#      that briefly existed here); and (b) being asked directly afterwards
#      how much of scipy was really needed, with a preference for not
#      depending on it at all if avoidable. Turned out to be exactly two
#      call sites (items 6 and 9 above), both scriptable by hand:
#        - _trilinear_interpolate(values_3d, origin, spacing, query_points):
#          a compact, standard 8-corner trilinear interpolation, written by
#          hand, replacing scipy.interpolate.RegularGridInterpolator for the
#          external-.cube case (build_potential_interpolator_from_cube()),
#          where the grid's (nx,ny,nz) layout is known for certain because
#          this codebase's own cube_reader.py controls the parsing. Verified
#          to reproduce scipy's RegularGridInterpolator output exactly
#          (bit-for-bit, in a side-by-side test against the synthetic test
#          .cube files from item 9) before scipy was removed.
#        - _nearest_neighbor_lookup(pts, vals, query_points, chunk_size):
#          a chunked, vectorised brute-force nearest-neighbour search,
#          replacing scipy.spatial.cKDTree for pDynamo3's OWN potential
#          grid (build_potential_interpolator()) specifically. Deliberately
#          NOT "upgraded" to reuse _trilinear_interpolate() by reshaping
#          pDynamo's flat gridValues into (nx,ny,nz): unlike the .cube case,
#          this codebase does not actually know, confirmed, whether such a
#          reshape would match pDynamo3's own internal storage order for
#          RegularGrid -- guessing wrong there would silently produce a
#          plausible-looking but WRONG interpolated surface, not a crash.
#          Nearest-neighbour works directly off the flat (point, value)
#          pairs with no ordering assumption at all, at some cost in
#          accuracy relative to true trilinear -- an intentional, explicit
#          trade-off in favour of correctness over precision here.
#      Net effect: EasyHybrid no longer depends on scipy for anything in
#      this file. The two SurfaceAnalysisWindow methods that use either
#      helper (build_potential_interpolator / build_potential_interpolator_
#      from_cube) are otherwise unchanged from item 6/9's description.
#
#  Caveat that applies to ALL of the pDynamo-facing code above (items 6, 9,
#  11): none of it could be executed against a real, running pDynamo3
#  installation in the environment this session's assistant had access to
#  (no pDynamo3 package was installed there). Everything that COULD be
#  exercised without pDynamo3 -- the .cube parser, the trilinear/nearest-
#  neighbour interpolators, the colormap math, the degenerate-triangle
#  filter -- was actually run and checked against synthetic data during the
#  session. The parts that call into pDynamo3's own compiled code
#  (RegularGrid, MarchingCubes_Isosurface3D, MakeVertexNormalsFromPolygonal
#  Normals, generator.Isosurface(), etc.) were instead validated by reading
#  pDynamo3's own source, cloned fresh from github.com/pdynamo/pDynamo3, and
#  were only actually confirmed working by the user running them live, in
#  their own environment where pDynamo3 is installed, over the course of
#  this same session.
#===============================================================================
#

import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk

#from GTKGUI.gtkWidgets.filechooser import FileChooser
#from easyhybrid.pDynamoMethods.pDynamo2Vismol import *
import gc
import os
import threading
import time
import numpy as np

from copy import deepcopy
import multiprocessing



from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from vismol.libgl.representations import SurfaceRepresentation
from vismol.core.vismol_object import VismolObject

#----------------------------------------------------------------------
from pSimulation import QCGridPropertyGenerator
from pCore       import *
from pScientific.Geometry3 import RegularGrid
from pScientific.Arrays     import Array, Reshape, RealArrayND
from pScientific.Surfaces   import MarchingCubes_Isosurface3D
# [EN] scipy (cKDTree / RegularGridInterpolator) was used here in an earlier
# version of this file for the MEP potential interpolation, but was DROPPED
# entirely: pDynamo3's own RegularGrid already exposes .origin / .spacing /
# .shape as plain Python properties (pScientific/Geometry3/__extensions__/
# pyrex/pScientific.Geometry3.RegularGrid.pyx), which is all that is needed
# to do genuine trilinear interpolation by hand, in pure numpy -- see
# _trilinear_interpolate() below, used by both build_potential_interpolator()
# (pDynamo grid) and build_potential_interpolator_from_cube() (external
# .cube grid). This removes a whole third-party dependency that had already
# caused a real problem: surface_analysis_window.py is imported
# UNCONDITIONALLY very early during EasyHybrid startup (from
# gui/main/main_window.py), and a machine without scipy installed got a
# ModuleNotFoundError all the way up through main_window.py -> gui.main ->
# easyhybrid.py, before the app even had a chance to open a window. As a
# bonus, this also replaces the pDynamo-grid nearest-neighbour lookup
# (cKDTree) with real trilinear interpolation -- strictly more accurate.
from util.colormaps import COLOR_MAPS
from util.cube_reader import read_cube_file, CubeFileError
#----------------------------------------------------------------------

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')




class SurfaceAnalysisWindow(Gtk.Window):
    """ Class doc """
    def __init__(self, main = None, system_liststore = None ):
        """ Class initialiser """
        self.main       = main#self.main.system_liststore
        self.home       = main.home
        self.visible    = False        
        self.p_session  = main.p_session
        self.vm_session = main.vm_session
        self.frame      = self.vm_session.frame
        
        self.orbital_liststore_dict = {
                                       #vm_object_id:[liststore1, liststore2...]
                                       }
        self.wave_function_dict = {
                                       #vm_object_id:[orbitals, gridgenerator]
                                      }
        self.counter = 1000
        
        self.modes = None
        self.running = False
        self.stop_thread = True
        
        self.lower = 0
        self.upper = 0
        
        self.selection_liststore       = Gtk.ListStore(str, int)
        self.selection_liststore_dict  =   {
                                           # system_e_id : Gtk.ListStore(str, int)
                                           }
    
    def open_window (self):
        """  """
        if self.visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/analysis/surface_analysis_window.glade')) #/home/fernando/programs/EasyHybrid3/gui/windows/normal_mode_analysis.glade
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_default_size(200, 600)  
            self.window.set_title('Surfaces')  
            self.window.set_keep_above(True)
            
            
            
            # - - - - - - - coordinates combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.coordinates_liststore = Gtk.ListStore(str, int, int)
            self.box2 = self.builder.get_object('box_coordinates')
            self.coordinates_combobox = CoordinatesComboBox(self.main.vobject_liststore_dict[self.p_session.active_id]) 
            self.box2.pack_start(self.coordinates_combobox, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''

            
            # - - - - - - - systems combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box1 = self.builder.get_object('box_system')
            self.system_names_combo = SystemComboBox (self.main, self.coordinates_combobox)
            #self.system_names_combo.connect("changed", self.on_system_names_combobox_changed)
            self.box1.pack_start(self.system_names_combo, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''

            

            self.btn_import_wfunction =  self.builder.get_object('btn_import_wavefunction')
            self.btn_import_wfunction.connect('clicked', self.on_button_import_wavefunction)


            #                       SURFACE TYPE COMBOBOX
            #'''--------------------------------------------------------------------------------------------'''
            self.box_surface_type = self.builder.get_object('box_surface_type')
            self.cbx_surface_type = Gtk.ComboBoxText()
            self.cbx_surface_type.connect("changed", self.surface_combobox_change)

            self.surface_options = ["Orbitals"               , 
                                    "Electrostatic Potential",
                                    "Density"              , 
                                    'External'              ,
                                    'MEP (density + potential)' ]
            for i, item in  enumerate(self.surface_options):
                self.cbx_surface_type.insert(i, str(i), item )
                
            self.box_surface_type.pack_start(self.cbx_surface_type, False, False, 0)
            #'''--------------------------------------------------------------------------------------------'''


            #                       RENDER MODE (superficie preenchida vs wireframe)
            #'''--------------------------------------------------------------------------------------------'''
            self.chk_surface_wireframe = Gtk.CheckButton(label="Wireframe")
            self.chk_surface_wireframe.connect("toggled", self.on_surface_wireframe_toggled)
            self.box_surface_type.pack_start(self.chk_surface_wireframe, False, False, 0)
            #'''--------------------------------------------------------------------------------------------'''


            #                       OPACIDADE (transparencia da superficie)
            #'''--------------------------------------------------------------------------------------------'''
            self.label_surface_opacity = Gtk.Label(label="Opacity:")
            self.scale_surface_opacity = Gtk.Scale.new_with_range ( Gtk.Orientation.HORIZONTAL, 0, 100, 1 )
            self.scale_surface_opacity.set_value ( 100 )   # 100% = opaco, igual ao comportamento de antes
            self.scale_surface_opacity.set_size_request ( 120, -1 )
            self.scale_surface_opacity.set_digits ( 0 )
            self.scale_surface_opacity.set_value_pos ( Gtk.PositionType.RIGHT )
            self.scale_surface_opacity.connect ( "value-changed", self.on_surface_opacity_changed )
            self.scale_surface_opacity.set_tooltip_text (
                "Opacidade das superficies (100% = opaco, 0% = totalmente "
                "transparente). Aplica em todas as superficies ja criadas "
                "nesta sessao (orbitais, potencial, densidade, MEP...).")

            self.box_surface_type.pack_start(self.label_surface_opacity, False, False, 0)
            self.box_surface_type.pack_start(self.scale_surface_opacity, True, True, 0)
            #'''--------------------------------------------------------------------------------------------'''


            #                       SHADING (flat vs smooth normals)
            #'''--------------------------------------------------------------------------------------------'''
            self.chk_surface_smooth = Gtk.CheckButton(label="Smooth shading")
            self.chk_surface_smooth.connect("toggled", self.on_surface_smooth_toggled)
            self.chk_surface_smooth.set_tooltip_text(
                "Desligado (padrao): normal constante por face (flat "
                "shading) -- mostra as facetas da triangulacao do marching "
                "cubes.\n"
                "Ligado: normal interpolada por vertice, media das faces "
                "adjacentes (smooth shading) -- superficie com aparencia "
                "mais lisa, sem facetas aparentes.")
            self.box_surface_type.pack_start(self.chk_surface_smooth, False, False, 0)
            #'''--------------------------------------------------------------------------------------------'''


            #                       MEP COLOR SCALE (vmin/vmax manuais, opcional)
            #'''--------------------------------------------------------------------------------------------'''
            self.label_mep_vmin = Gtk.Label(label="MEP vmin:")
            self.entry_mep_vmin = Gtk.Entry()
            self.entry_mep_vmin.set_width_chars(8)
            self.entry_mep_vmin.set_placeholder_text("auto")
            self.entry_mep_vmin.set_tooltip_text(
                "Valor minimo do potencial eletrostatico (em unidades "
                "atomicas, hartree/e) usado na escala de cor do MEP.\n"
                "Vazio = automatico (percentil 2% dos valores calculados).")

            self.label_mep_vmax = Gtk.Label(label="vmax:")
            self.entry_mep_vmax = Gtk.Entry()
            self.entry_mep_vmax.set_width_chars(8)
            self.entry_mep_vmax.set_placeholder_text("auto")
            self.entry_mep_vmax.set_tooltip_text(
                "Valor maximo do potencial eletrostatico (em unidades "
                "atomicas, hartree/e) usado na escala de cor do MEP.\n"
                "Vazio = automatico (percentil 98% dos valores calculados).")

            self.box_surface_type.pack_start(self.label_mep_vmin, False, False, 0)
            self.box_surface_type.pack_start(self.entry_mep_vmin, False, False, 0)
            self.box_surface_type.pack_start(self.label_mep_vmax, False, False, 0)
            self.box_surface_type.pack_start(self.entry_mep_vmax, False, False, 0)
            # comecam escondidos -- so fazem sentido com "MEP" selecionado
            # (index == 4), ver surface_combobox_change().
            self.label_mep_vmin.hide()
            self.entry_mep_vmin.hide()
            self.label_mep_vmax.hide()
            self.entry_mep_vmax.hide()
            #'''--------------------------------------------------------------------------------------------'''


            #                       MEP COLORMAP (COLOR_MAPS de util/colormaps.py)
            #'''--------------------------------------------------------------------------------------------'''
            self.label_mep_cmap = Gtk.Label(label="Colormap:")
            self.cbx_mep_cmap   = Gtk.ComboBoxText()
            self.cbx_mep_cmap.set_tooltip_text(
                "Colormap usado pra mapear o potencial eletrostatico em cor "
                "(ver COLOR_MAPS em util/colormaps.py). Mapas divergentes "
                "(coolwarm, vik, berlin, bam) tem convencao vermelho="
                "negativo/azul=positivo; mapas sequenciais (jet, rainbow, "
                "gnuplot, magma, viridis...) nao tem centro definido.")

            self._mep_cmap_names = list ( COLOR_MAPS.keys() )
            for i, name in enumerate ( self._mep_cmap_names ):
                self.cbx_mep_cmap.insert ( i, str(i), name )
            try:
                default_idx = self._mep_cmap_names.index ( 'coolwarm' )
            except ValueError:
                default_idx = 0
            self.cbx_mep_cmap.set_active ( default_idx )

            self.box_surface_type.pack_start(self.label_mep_cmap, False, False, 0)
            self.box_surface_type.pack_start(self.cbx_mep_cmap, False, False, 0)
            self.label_mep_cmap.hide()
            self.cbx_mep_cmap.hide()
            #'''--------------------------------------------------------------------------------------------'''


            self.label_frame = self.builder.get_object('label_frame')
                 
            system  = self.main.p_session.get_system()

            for row in self.main.vobject_liststore_dict[self.main.p_session.active_id]:
                a = list (row)
                self.coordinates_liststore.append([a[0], a[1], a[2]])

            self.coordinates_combobox.set_model(self.coordinates_liststore)
            

            columns = [' ', 'Orbital', 'Occ.', 'Energy']#, 'visible']
            
            self.liststore = Gtk.ListStore(int, str, int, float)#, bool)

            self.treeview = self.builder.get_object('selection_treeview')#Gtk.TreeView(model=self.liststore)
            
            # Remove todas as colunas criadas pelo Glade
            for col in self.treeview.get_columns():
                print(col)
                self.treeview.remove_column(col)
            
            self.treeview.set_model(self.liststore)
             
            for i, column_title in enumerate(columns):
                renderer = Gtk.CellRendererText()
                #if column_title == 'visible':
                #    renderer_toggle = Gtk.CellRendererToggle()
                #    renderer_toggle.connect("toggled", self.on_cell_toggled)
                #    column = Gtk.TreeViewColumn(column_title, renderer_toggle, active=4)
                #else:    
                column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                self.treeview.append_column(column)
            
            #-----------------------------------------------------------
            self.btn_render = self.builder.get_object('btn_render')
            self.btn_render.connect('clicked', self.on_render_button)
            #-----------------------------------------------------------
            
            
            self.btn_color_plus  = self.builder.get_object('btn_color_plus')
            self.btn_color_minus = self.builder.get_object('btn_color_minus')
            # Definindo uma cor específica (por exemplo, vermelho)
            rgba = Gdk.RGBA()
            
            rgba.parse("blue")  # ou "rgb(255,0,0)", ou "#FF0000"
            self.btn_color_minus.set_rgba(rgba)
            
            rgba.parse("red") 
            self.btn_color_plus.set_rgba(rgba)
            
            
            try:
                active_e_id = self.p_session.active_id
                self.system_names_combo.set_active_system(e_id = active_e_id) 
            except:
                pass
            
            #self.refresh_system_liststore()
            #self.treeview_menu         = TreeViewMenu(self)
            
            self.window.show_all()                                               
            
            
            #                       EXTERNAL CUBE FILES (density/orbital obrigatorio, potential opcional p/ MEP)
            #'''--------------------------------------------------------------------------------------------'''
            self.external_density_path   = None
            self.external_potential_path = None

            self.btn_external_file = self.builder.get_object('btn_external_file')
            self.btn_external_file.connect("file-set", self.on_external_density_file_set)
            self.btn_external_file.set_tooltip_text(
                "Arquivo .cube de densidade ou orbital (ex: gerado pelo orca_plot "
                "do ORCA). A malha da superficie vem daqui.")

            self.label_external_potential = Gtk.Label(label="Potential .cube (opcional p/ MEP):")
            self.btn_external_potential_file = Gtk.FileChooserButton(
                title  = "Selecione o arquivo .cube de potencial (opcional)",
                action = Gtk.FileChooserAction.OPEN,
            )
            _cube_filter = Gtk.FileFilter()
            _cube_filter.set_name("Arquivos Cube (*.cube)")
            _cube_filter.add_pattern("*.cube")
            self.btn_external_potential_file.add_filter(_cube_filter)
            self.btn_external_potential_file.connect("file-set", self.on_external_potential_file_set)
            self.btn_external_potential_file.set_tooltip_text(
                "Opcional. Se fornecido, colore a malha de densidade/orbital "
                "acima por potencial eletrostatico interpolado desse cubo "
                "(mapa continuo de cor, igual ao MEP -- ver mep_colormap).")
            self.box_surface_type.pack_start(self.label_external_potential, False, False, 0)
            self.box_surface_type.pack_start(self.btn_external_potential_file, False, False, 0)
            self.label_external_potential.hide()
            self.btn_external_potential_file.hide()
            #'''--------------------------------------------------------------------------------------------'''

            # so agora, com TODOS os widgets acima ja criados, e que ativamos o
            # item 0 do combobox -- set_active() dispara "changed" (chama
            # surface_combobox_change) IMEDIATAMENTE e de forma sincrona, entao
            # precisa vir depois de tudo que esse handler pode tentar
            # mostrar/esconder (senao da AttributeError: widget ainda nao existe).
            self.cbx_surface_type.set_active(0)

            self.builder.get_object('btn_external_file').hide()
            self.builder.get_object('label_external_file').hide()
            #self.system_names_combo.set_active(0)
            self.visible    =  True
            '''--------------------------------------------------------------------------------------------'''

    def close_window (self, button, data  = None):
        """ Function doc """

        #self.stop(None)
        #self.running = False  


        self.window.destroy()
        self.visible    =  False
        #print('self.visible',self.visible)

    def on_cell_toggled(self, widget, path):
        self.liststore[path][4] = not self.liststore[path][4]

    def update_window (self, system_names = True, coordinates = False,  selections = True):
        """ Function doc """
        if self.visible:
            
            _id = self.system_names_combo.get_active()
            if _id == -1:
                '''_id = -1 means no item inside the combobox'''
                #self.selection_liststore.clear()
                #self.coordinates_liststore.clear()
                return None
            else:    
                _, system_id = self.main.system_liststore[_id]
            
            
            #if system_names:
            #    self.refresh_system_liststore ()
            #    self.system_names_combo.set_active(_id)
            #
            #if coordinates:
            #    self.refresh_coordinates_liststore ()
            #
            #
            if selections:
                _, system_id = self.main.system_liststore[_id]
                self.refresh_selection_liststore(system_id)
        else:
            pass
     
    def update (self, system_names = True, coordinates = False,  selections = True ):
        """ Function doc """
        pass
    
    def refresh_system_liststore (self):
        """ Function doc """
        self.main.refresh_system_liststore()

    # [EN] The next few handlers (wireframe / opacity / smooth-shading) all
    # follow the same pattern: walk every vismol object flagged
    # is_surface == True in self.vm_session.vm_objects_dic, and call one
    # setter method on each of its representations. None of them regenerate
    # the mesh -- they only flip a rendering-time flag/uniform that was
    # added to SurfaceRepresentation (representations.py) during this
    # session. See changelog items 2, 3 and 4 at the top of this file.
    def on_surface_wireframe_toggled (self, widget):
        """ Alterna entre superficie preenchida (GL_FILL) e wireframe (GL_LINE)
        para todas as representacoes de superficie ja criadas nesta sessao.
        Nao recalcula malha nem recompila shader -- so muda o render_mode
        de cada SurfaceRepresentation encontrada (ver representations.py). """
        mode = "lines" if widget.get_active() else "surface"
        for vobject in self.vm_session.vm_objects_dic.values():
            if not getattr(vobject, "is_surface", False):
                continue
            for rep in vobject.representations.values():
                if hasattr(rep, "set_render_mode"):
                    rep.set_render_mode(mode)
        self.vm_session.vm_glcore.queue_draw()

    def on_surface_opacity_changed (self, widget):
        """ Ajusta a opacidade (alpha) de todas as representacoes de
        superficie ja criadas nesta sessao. widget.get_value() vai de 0
        a 100 (%); SurfaceRepresentation.set_alpha() espera 0.0-1.0. """
        alpha = widget.get_value() / 100.0
        for vobject in self.vm_session.vm_objects_dic.values():
            if not getattr(vobject, "is_surface", False):
                continue
            for rep in vobject.representations.values():
                if hasattr(rep, "set_alpha"):
                    rep.set_alpha(alpha)
        self.vm_session.vm_glcore.queue_draw()

    def on_surface_smooth_toggled (self, widget):
        """ Alterna entre flat shading (normal por face) e smooth shading
        (normal por vertice, media das faces adjacentes) em todas as
        representacoes de superficie ja criadas nesta sessao. Nao
        recalcula a malha -- as normais suaves ja foram calculadas na
        geracao (surface.MakeVertexNormalsFromPolygonalNormals(), nativo
        do pDynamo3) e enviadas pro VAO; aqui so muda qual delas o
        shader usa (ver geometry_shader_surface). """
        mode = "smooth" if widget.get_active() else "flat"
        for vobject in self.vm_session.vm_objects_dic.values():
            if not getattr(vobject, "is_surface", False):
                continue
            for rep in vobject.representations.values():
                if hasattr(rep, "set_shading_mode"):
                    rep.set_shading_mode(mode)
        self.vm_session.vm_glcore.queue_draw()

    # [EN] btn_external_file is a GtkFileChooserButton defined in the .glade
    # file (NOT a plain Gtk.Button) -- it only emits "file-set", never
    # "clicked". An earlier attempt to connect("clicked", ...) on it crashed
    # live with "TypeError: unknown signal name: clicked" the first time
    # this code actually ran outside this environment. See changelog item 9
    # / 10(b). btn_external_potential_file is a second FileChooserButton,
    # this one created in pure Python for symmetry, for the optional
    # potential cube used to MEP-colour the imported surface.
    def on_external_density_file_set (self, widget):
        """ Chamado quando o usuario escolhe o .cube de densidade/orbital
        no proprio GtkFileChooserButton (widget nativo do .glade --
        emite 'file-set', nao 'clicked'; um botao comum nao serve aqui). """
        path = widget.get_filename ( )
        if path is not None:
            self.external_density_path = path

    def on_external_potential_file_set (self, widget):
        """ Chamado quando o usuario escolhe o .cube de potencial
        (opcional -- se fornecido, colore a malha por MEP em vez de cor
        uniforme). """
        path = widget.get_filename ( )
        if path is not None:
            self.external_potential_path = path

    def surface_combobox_change (self, widget):
        """ Function doc """
        index = self.cbx_surface_type.get_active()
        print(index)

        if index == 4:
            self.label_mep_vmin.show()
            self.entry_mep_vmin.show()
            self.label_mep_vmax.show()
            self.entry_mep_vmax.show()
            self.label_mep_cmap.show()
            self.cbx_mep_cmap.show()
        else:
            self.label_mep_vmin.hide()
            self.entry_mep_vmin.hide()
            self.label_mep_vmax.hide()
            self.entry_mep_vmax.hide()
            self.label_mep_cmap.hide()
            self.cbx_mep_cmap.hide()

        if index in [1,2,4]:
            #self.builder.get_object('btn_import_wavefunction').set_sensitive(False)
            #self.builder.get_object('selection_treeview')     .set_sensitive(False)
            self.builder.get_object('btn_import_wavefunction').hide()
            self.builder.get_object('selection_treeview')     .hide()
            self.builder.get_object('label_external_file').hide()
            self.builder.get_object('btn_external_file')  .hide()
            self.label_external_potential.hide()
            self.btn_external_potential_file.hide()
        
        elif index == 3:
            self.builder.get_object('label_external_file').show()
            self.builder.get_object('btn_external_file').show()
            self.label_external_potential.show()
            self.btn_external_potential_file.show()
            self.builder.get_object('btn_import_wavefunction').set_sensitive(False)
            self.builder.get_object('selection_treeview'     ).set_sensitive(False)
            self.builder.get_object('btn_import_wavefunction').hide()
            self.builder.get_object('selection_treeview'     ).hide()
        
        else:
            self.builder.get_object('btn_import_wavefunction').set_sensitive(True)
            self.builder.get_object('selection_treeview')     .set_sensitive(True)
            self.builder.get_object('btn_import_wavefunction').show()
            self.builder.get_object('selection_treeview')     .show()
            self.builder.get_object('label_external_file').hide()
            self.builder.get_object('btn_external_file')  .hide()
            self.label_external_potential.hide()
            self.btn_external_potential_file.hide()

    def on_coordinates_combobox_changed(self, widget):
        """ Function doc """
        index = self.coordinates_combobox.get_active()
       
        if index == -1:
            '''_id = -1 means no item inside the combobox'''
            return None
        
        else:    
            name, vobject_id, system_id = self.coordinates_liststore[index]
            
            #print(name, vobject_id, system_id)
            
            try:
                vobject = self.vm_session.vm_objects_dic[vobject_id]
                self.liststore.clear()
                for i, data in vobject.normal_modes_dict.items():
                    self.liststore.append([False, str(i), data[0]])
            except:
                print('vobject has no Normal Modes data')
                pass

    def on_treeview_Objects_row_activated(self, tree, event, data):
        rgba_plus = self.btn_color_minus.get_rgba()
        rgba_minus = self.btn_color_plus.get_rgba()
        color_plus  = [rgba_plus.red , rgba_plus.green,  rgba_plus.blue ]
        color_minus = [rgba_minus.red, rgba_minus.green, rgba_minus.blue]
        
        system_id = self.system_names_combo.get_system_id()
        system    = self.main.p_session.psystem[system_id]

        vobject_id    = self.coordinates_combobox.get_vobject_id()
        vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]

        _isovalue    = float(self.builder.get_object('entry_isovalue').get_text())
        _GridSpacing = float(self.builder.get_object('entry_spacing') .get_text())
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        

        #system  = self.main.p_session.get_system()
        #'''
        backup = []
        try:
            backup.append(system.e_treeview_iter)
            backup.append(system.e_liststore_iter)
            system.e_treeview_iter   = None
            system.e_liststore_iter  = None
        except:
            pass



        '''
        "key" is the acesses key to the dictionary containg the selection lists
        there is no two selection lists with the same name.
        indexes =  A list of atoms for selection
        '''
        key     = model.get_value(iter, 0)
        print(key, vismol_object.frames.shape[0], model, model[iter][1])
        name = str(key) +' '+model[iter][1]#+' '+ str(model[iter][3])
        #_GridSpacing = 0.6
        _OrbitalTag    = "Grid Orbitals"
        _IsosurfaceTag = "Isosurface"
        
        
        trajectory = [None]*vismol_object.frames.shape[0]
        joblist = []
        
        for frame in range(vismol_object.frames.shape[0]):
            #self.p_session.set_psystem_coordinates_from_vobject(vobject)
            #'''
            print(vismol_object, frame)
            self.p_session.set_psystem_coordinates_from_vobject( vobject = vismol_object, 
                                                                           system_id = None, 
                                                                           frame = frame)
            
            parameters = {
            'type'           : 'orbital',
            '_GridSpacing'   : _GridSpacing,
            '_OrbitalTag'    : _OrbitalTag,
            '_isovalue'      : _isovalue,
            '_IsosurfaceTag' : _IsosurfaceTag,
            'orbital_key'    : key,
            'color_plus'     : color_plus ,
            'color_minus'    : color_minus,
            }

            coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
            
            joblist.append([frame, system, coords, parameters])
        
        p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
        results = p.map(generate_grid_parallel, joblist)
        
        #if interface:
        try:
            system.e_treeview_iter   = backup[0]
            system.e_liststore_iter  = backup[1]
        except:
            pass
        
        #vobject_tmp = VismolObject(name="UNK", index=-1,
        vobject_tmp = VismolObject(name= name , index=-1,
                                   vismol_session        = self.vm_session,
                                   trajectory            = [],
                                   bonds_pair_of_indexes = [0,1])
        
        
        vobject_tmp.model_mat = vismol_object.model_mat 
        vobject_tmp.trans_mat = vismol_object.trans_mat 
        
        #vismol_object.surface_trajectory = results # trajectory
        vobject_tmp.surface_trajectory = results # trajectory
        #-----------------------------------------------------------------------
        #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
        #vismol_object.representations["surface1"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
        vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                           vismol_glcore = self.vm_session.vm_glcore ,  
                                                                           name          = 'surface'                 ,
                                                                           active        = True                      ,
                                                                           indexes       = []                        ,
                                                                           is_dynamic    = False                     ,
                                                                           iso_color     = [1,0,0]                   ,
                                                                           surface_name  = 'obital_plus'                )
                                                     

        #-----------------------------------------------------------------------
        #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
        #vismol_object.representations["surface2"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
        vobject_tmp.representations["surface2"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                           vismol_glcore = self.vm_session.vm_glcore ,  
                                                                           name          = 'surface'                 ,
                                                                           active        = True                      ,
                                                                           indexes       = []                        ,
                                                                           is_dynamic    = False                     ,
                                                                           iso_color     = [0,0,1]                   ,
                                                                           surface_name  = 'obital_minus'           )
        
        vobject_tmp.parameters = parameters
        
        vobject_tmp.frames = vismol_object.frames
        vobject_tmp.active = True
        vobject_tmp.is_surface = True
        vobject_tmp.e_id = system.e_id
        self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
        
        
        self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
        self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
        self.main.refresh_widgets()
        self.vm_session.vm_glcore.queue_draw()
        self.counter +=1
        
        
        
        '''
        for frame , data in enumerate(self.wave_function_dict[vobject_id]):
            
            
            generator = data[2]
            #system.Energy()
            parameters = {
            '_GridSpacing'   : _GridSpacing,
            '_OrbitalTag'    : _OrbitalTag,
            '_isovalue'      : _isovalue,
            '_IsosurfaceTag' : _IsosurfaceTag,
            'orbital_key'    : key,
            }
            
            joblist.append([frame, generator, parameters])
        
        p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
        #'''
        '''
        results = p.map(generate_grid_parallel, joblist)
        
        print (results)
        vismol_object.surface_trajectory = results # trajectory
        #-----------------------------------------------------------------------
        #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
        vismol_object.representations["surface1"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
                                                                           vismol_glcore = self.vm_session.vm_glcore ,  
                                                                           name          = 'surface'                 ,
                                                                           active        = True                      ,
                                                                           indexes       = []                        ,
                                                                           is_dynamic    = False                     ,
                                                                           iso_color     = [1,0,0]                   ,
                                                                           surface_name  = 'obital_plus'                )
                                                     

        #-----------------------------------------------------------------------
        #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
        vismol_object.representations["surface2"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
                                                                           vismol_glcore = self.vm_session.vm_glcore ,  
                                                                           name          = 'surface'                 ,
                                                                           active        = True                      ,
                                                                           indexes       = []                        ,
                                                                           is_dynamic    = False                     ,
                                                                           iso_color     = [0,0,1]                   ,
                                                                           surface_name  = 'obital_minus'           )
        self.vm_session.vm_glcore.queue_draw()
        #'''
    
    def _update_liststore (self):
        """ Function doc """
        
        vobject_id = self.coordinates_combobox.get_vobject_id()
        print(vobject_id)
        
        model = self.treeview.get_model()
        if model is not None:
            # Remove todos os itens do modelo
            model.clear()
       
        if self.frame > len(self.orbital_liststore_dict[vobject_id]):
            self.treeview.set_model(self.orbital_liststore_dict[vobject_id][-1])
        else:
            
            #for frame , data in enumerate(self.wave_function_dict[vobject_id]):
            orbitals = self.wave_function_dict[vobject_id][self.frame][0]
            for i in range(len(orbitals)):
                reverse_index = -i-1 #- len(orbitals)
                model.append(orbitals[reverse_index ])

    def set_frame (self ):
        """ Function doc """
        self.frame =  self.vm_session.frame
        self.label_frame.set_text('Frame = {}'.format(self.frame))
        self._update_liststore()

    def on_render_button (self, widget):
 
 
        rgba_plus = self.btn_color_minus.get_rgba()
        rgba_minus = self.btn_color_plus.get_rgba()
        color_plus  = [rgba_plus.red , rgba_plus.green,  rgba_plus.blue ]
        color_minus = [rgba_minus.red, rgba_minus.green, rgba_minus.blue]
        
        index = self.cbx_surface_type.get_active()
        print(index, color_minus, color_plus)

        system_id = self.system_names_combo.get_system_id()
        system    = self.main.p_session.psystem[system_id]

        vobject_id    = self.coordinates_combobox.get_vobject_id()
        vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]

        _isovalue    = float(self.builder.get_object('entry_isovalue').get_text())
        _GridSpacing = float(self.builder.get_object('entry_spacing') .get_text())
        
        
        backup = []
        try:
            backup.append(system.e_treeview_iter)
            backup.append(system.e_liststore_iter)
            system.e_treeview_iter   = None
            system.e_liststore_iter  = None
        except:
            pass
        
        
        if index == 2:
            joblist = []
            for frame in range(vismol_object.frames.shape[0]):
                #'''
                self.p_session.set_psystem_coordinates_from_vobject( vobject = vismol_object, 
                                                                               system_id = None, 
                                                                               frame = frame)
                parameters = {
                'type'           : 'density',
                '_GridSpacing'   : _GridSpacing,
                '_OrbitalTag'    : 'density',
                '_isovalue'      : _isovalue,
                '_IsosurfaceTag' : 'density',
                'orbital_key'    : 0,
                'color_plus'     : color_plus  ,
                'color_minus'    : color_minus ,
                
                }
                coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
                joblist.append([frame, system, coords, parameters])
                
            p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
            results = p.map(generate_grid_parallel, joblist)
            
            #if interface:
            try:
                system.e_treeview_iter   = backup[0]
                system.e_liststore_iter  = backup[1]
            except:
                pass
            name ='Density'
            #vobject_tmp = VismolObject(name="UNK", index=-1,
            vobject_tmp = VismolObject(name= name , index=-1,
                                        vismol_session        = self.vm_session,
                                        trajectory            = [],
                                        bonds_pair_of_indexes = [0,1])
            
            
            vobject_tmp.model_mat = vismol_object.model_mat 
            vobject_tmp.trans_mat = vismol_object.trans_mat 
            
            #vismol_object.surface_trajectory = results # trajectory
            vobject_tmp.surface_trajectory = results # trajectory
            #-----------------------------------------------------------------------
            #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
            #vismol_object.representations["surface1"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
            vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                                vismol_glcore = self.vm_session.vm_glcore ,  
                                                                                name          = 'surface'                 ,
                                                                                active        = True                      ,
                                                                                indexes       = []                        ,
                                                                                is_dynamic    = False                     ,
                                                                                iso_color     = color_plus                   ,
                                                                                surface_name  = 'obital_plus'                )
                                                            
        
            #-----------------------------------------------------------------------
            #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
            #vismol_object.representations["surface2"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
            
            #vobject_tmp.representations["surface2"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
            #                                                                    vismol_glcore = self.vm_session.vm_glcore ,  
            #                                                                    name          = 'surface'                 ,
            #                                                                    active        = True                      ,
            #                                                                    indexes       = []                        ,
            #                                                                    is_dynamic    = False                     ,
            #                                                                    iso_color     = [0,0,1]                   ,
            #                                                                    surface_name  = 'obital_minus'           )
            #
            
            vobject_tmp.parameters = parameters
            
            vobject_tmp.frames = vismol_object.frames
            vobject_tmp.active = False
            vobject_tmp.is_surface = True
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            print('\n\nvismol_object.e_treeview_iter', vismol_object.e_treeview_iter,'\n\n')
            self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
            # Add the VisMol object to the vobject liststore dictionary
            self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
            # Refresh the widgets in the main window
            self.main.refresh_widgets()
            
            
            
            
            #print(vobject_tmp, vobject_tmp.surface_trajectory)
            #self.vm_session.vm_objects_dic[self.counter] = vobject_tmp
            self.vm_session.vm_glcore.queue_draw()
            self.counter +=1

        elif index == 4:
            def _parse_optional_float ( entry ):
                text = entry.get_text().strip()
                if text == "":
                    return None
                try:
                    return float ( text )
                except ValueError:
                    return None   # texto invalido -- cai no automatico (percentil)

            _mep_vmin = _parse_optional_float ( self.entry_mep_vmin )
            _mep_vmax = _parse_optional_float ( self.entry_mep_vmax )
            _mep_cmap_idx = self.cbx_mep_cmap.get_active()
            _mep_cmap_name = ( self._mep_cmap_names[_mep_cmap_idx]
                                if _mep_cmap_idx >= 0 else 'coolwarm' )

            joblist = []
            for frame in range(vismol_object.frames.shape[0]):
                self.p_session.set_psystem_coordinates_from_vobject( vobject = vismol_object, 
                                                                               system_id = None, 
                                                                               frame = frame)
                parameters = {
                'type'           : 'mep',
                '_GridSpacing'   : _GridSpacing,
                '_OrbitalTag'    : 'density_mep',
                '_isovalue'      : _isovalue,
                '_IsosurfaceTag' : 'Isosurface',
                'orbital_key'    : 0,
                'color_plus'     : color_plus  ,
                'color_minus'    : color_minus ,
                'mep_vmin'       : _mep_vmin   ,
                'mep_vmax'       : _mep_vmax   ,
                'mep_cmap_name'  : _mep_cmap_name ,
                }
                coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
                joblist.append([frame, system, coords, parameters])
                
            p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
            results = p.map(generate_grid_parallel, joblist)
            
            try:
                system.e_treeview_iter   = backup[0]
                system.e_liststore_iter  = backup[1]
            except:
                pass
            name ='MEP'
            vobject_tmp = VismolObject(name= name , index=-1,
                                        vismol_session        = self.vm_session,
                                        trajectory            = [],
                                        bonds_pair_of_indexes = [0,1])
            
            vobject_tmp.model_mat = vismol_object.model_mat 
            vobject_tmp.trans_mat = vismol_object.trans_mat 
            vobject_tmp.surface_trajectory = results # trajectory
            #-----------------------------------------------------------------------
            # iso_color aqui e so um placeholder: a cor de verdade ja vem
            # por vertice dentro de surface_trajectory (calculada pelo
            # mep_colormap em generate_grid_parallel), nao usada por
            # SurfaceRepresentation.__init__ (self.iso_color nunca e
            # armazenado -- ver representations.py).
            vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                                vismol_glcore = self.vm_session.vm_glcore ,  
                                                                                name          = 'surface'                 ,
                                                                                active        = True                      ,
                                                                                indexes       = []                        ,
                                                                                is_dynamic    = False                     ,
                                                                                iso_color     = color_plus                ,
                                                                                surface_name  = 'obital_plus'                )
            
            vobject_tmp.parameters = parameters
            
            vobject_tmp.frames = vismol_object.frames
            vobject_tmp.active = False
            vobject_tmp.is_surface = True
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
            self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
            self.main.refresh_widgets()
            self.vm_session.vm_glcore.queue_draw()
            self.counter +=1

        elif index == 3:
            if not self.external_density_path:
                print("Nenhum arquivo .cube de densidade/orbital selecionado -- "
                      "clique em 'Escolher arquivo...' antes de renderizar.")
                return False

            _mep_cmap_idx = self.cbx_mep_cmap.get_active() if hasattr(self, "_mep_cmap_names") else -1
            _mep_cmap_name = ( self._mep_cmap_names[_mep_cmap_idx]
                                if _mep_cmap_idx >= 0 else 'coolwarm' )

            parameters = {
            'type'                     : 'external_cube',
            '_isovalue'                : _isovalue,
            'color_plus'               : color_plus  ,
            'color_minus'              : color_minus ,
            'external_density_path'   : self.external_density_path   ,
            'external_potential_path' : self.external_potential_path , # None = sem MEP, cor uniforme
            'mep_vmin'      : None,   # sem campo manual pra External ainda -- sempre automatico (percentil)
            'mep_vmax'      : None,
            'mep_cmap_name' : _mep_cmap_name,
            }

            # cubo externo e um arquivo estatico -- sem system/coords/QC
            # nenhum envolvido, so leitura de arquivo + marching cubes.
            # Roda direto (sem multiprocessing.Pool: e so I/O + um
            # algoritmo compilado, nao ha calculo QC pesado a paralelizar
            # aqui, e cada mudanca de aba abriria um Pool novo a toa).
            try:
                single_result = generate_grid_parallel ( [ 0, None, None, parameters ] )
            except CubeFileError as error:
                print ( "Erro lendo arquivo .cube: {}".format ( error ) )
                return False

            # replica o mesmo resultado pra todos os "frames" do objeto
            # pai, so por seguranca (surface_trajectory[frame] nao pode
            # dar index error se o usuario tiver uma trajetoria carregada
            # e trocar de frame -- o cubo externo e sempre a mesma
            # malha estatica, nao muda por frame).
            results = [ single_result ] * max ( 1, vismol_object.frames.shape[0] )

            try:
                system.e_treeview_iter   = backup[0]
                system.e_liststore_iter  = backup[1]
            except:
                pass
            name = 'External Cube'
            vobject_tmp = VismolObject(name= name , index=-1,
                                        vismol_session        = self.vm_session,
                                        trajectory            = [],
                                        bonds_pair_of_indexes = [0,1])
            
            vobject_tmp.model_mat = vismol_object.model_mat 
            vobject_tmp.trans_mat = vismol_object.trans_mat 
            vobject_tmp.surface_trajectory = results # trajectory
            vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                                vismol_glcore = self.vm_session.vm_glcore ,  
                                                                                name          = 'surface'                 ,
                                                                                active        = True                      ,
                                                                                indexes       = []                        ,
                                                                                is_dynamic    = False                     ,
                                                                                iso_color     = color_plus                ,
                                                                                surface_name  = 'obital_plus'                )
            
            vobject_tmp.parameters = parameters
            
            vobject_tmp.frames = vismol_object.frames
            vobject_tmp.active = False
            vobject_tmp.is_surface = True
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
            self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
            self.main.refresh_widgets()
            self.vm_session.vm_glcore.queue_draw()
            self.counter +=1

        elif index ==1:
            joblist = []
            for frame in range(vismol_object.frames.shape[0]):
                #'''
                self.p_session.set_psystem_coordinates_from_vobject( vobject = vismol_object, 
                                                                               system_id = None, 
                                                                               frame = frame)
                parameters = {
                'type'           : 'potential',
                '_GridSpacing'   : _GridSpacing,
                '_OrbitalTag'    : 'potential',
                '_isovalue'      : _isovalue,
                '_IsosurfaceTag' : 'Isosurface',
                'orbital_key'    : 0,
                'color_plus'     : color_plus  ,
                'color_minus'    : color_minus ,
                }
                coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
                joblist.append([frame, system, coords, parameters])
                
            p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
            results = p.map(generate_grid_parallel, joblist)
            
            #if interface:
            try:
                system.e_treeview_iter   = backup[0]
                system.e_liststore_iter  = backup[1]
            except:
                pass
            name ='Potential'
            vobject_tmp = VismolObject(name= name , index=-1,
                                        vismol_session        = self.vm_session,
                                        trajectory            = [],
                                        bonds_pair_of_indexes = [0,1])
            
            
            vobject_tmp.model_mat = vismol_object.model_mat 
            vobject_tmp.trans_mat = vismol_object.trans_mat 
            vobject_tmp.surface_trajectory = results # trajectory
            #-----------------------------------------------------------------------
            vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                                vismol_glcore = self.vm_session.vm_glcore ,  
                                                                                name          = 'surface'                 ,
                                                                                active        = True                      ,
                                                                                indexes       = []                        ,
                                                                                is_dynamic    = False                     ,
                                                                                iso_color     = color_plus                 ,
                                                                                surface_name  = 'obital_plus'             )
                                                            
        
            #-----------------------------------------------------------------------
            vobject_tmp.representations["surface2"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                                vismol_glcore = self.vm_session.vm_glcore ,  
                                                                                name          = 'surface'                 ,
                                                                                active        = True                      ,
                                                                                indexes       = []                        ,
                                                                                is_dynamic    = False                     ,
                                                                                iso_color     = color_minus               ,
                                                                                surface_name  = 'obital_minus'           )
            
            
            vobject_tmp.parameters = parameters
            
            vobject_tmp.frames = vismol_object.frames
            vobject_tmp.active = False
            vobject_tmp.is_surface = True
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            
            self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
            # Add the VisMol object to the vobject liststore dictionary
            self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
            # Refresh the widgets in the main window
            self.main.refresh_widgets()
            self.vm_session.vm_glcore.queue_draw()
            self.counter +=1

        elif index == 0:
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
        
            
            backup = []
            try:
                backup.append(system.e_treeview_iter)
                backup.append(system.e_liststore_iter)
                system.e_treeview_iter   = None
                system.e_liststore_iter  = None
            except:
                pass



            '''
            "key" is the acesses key to the dictionary containg the selection lists
            there is no two selection lists with the same name.
            indexes =  A list of atoms for selection
            '''
            key     = model.get_value(iter, 0)
            print(key, vismol_object.frames.shape[0], model, model[iter][1])
            name = str(key) +' '+model[iter][1]#+' '+ str(model[iter][3])
            #_GridSpacing = 0.6
            _OrbitalTag    = "Grid Orbitals"
            _IsosurfaceTag = "Isosurface"
            
            
            trajectory = [None]*vismol_object.frames.shape[0]
            joblist = []
            
            for frame in range(vismol_object.frames.shape[0]):
                #self.p_session.set_psystem_coordinates_from_vobject(vobject)
                #'''
                self.p_session.set_psystem_coordinates_from_vobject( vobject = vismol_object, 
                                                                               system_id = None, 
                                                                               frame = frame)
                
                parameters = {
                'type'           : 'orbital',
                '_GridSpacing'   : _GridSpacing,
                '_OrbitalTag'    : _OrbitalTag,
                '_isovalue'      : _isovalue,
                '_IsosurfaceTag' : _IsosurfaceTag,
                'orbital_key'    : key,
                'color_plus'     : color_plus  ,
                'color_minus'    : color_minus ,
                }

                coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
                
                joblist.append([frame, system, coords, parameters])
            
            p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
            results = p.map(generate_grid_parallel, joblist)
            
            #if interface:
            try:
                system.e_treeview_iter   = backup[0]
                system.e_liststore_iter  = backup[1]
            except:
                pass
            
            vobject_tmp = VismolObject(name= name , index=-1,
                                       vismol_session        = self.vm_session,
                                       trajectory            = [],
                                       bonds_pair_of_indexes = [0,1])
            
            vobject_tmp.model_mat = vismol_object.model_mat 
            vobject_tmp.trans_mat = vismol_object.trans_mat 
            vobject_tmp.surface_trajectory = results # trajectory
            #-----------------------------------------------------------------------
            vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                               vismol_glcore = self.vm_session.vm_glcore ,  
                                                                               name          = 'surface'                 ,
                                                                               active        = True                      ,
                                                                               indexes       = []                        ,
                                                                               is_dynamic    = False                     ,
                                                                               iso_color     = color_plus                   ,
                                                                               surface_name  = 'obital_plus'                )
                                                         

            #-----------------------------------------------------------------------
            vobject_tmp.representations["surface2"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                               vismol_glcore = self.vm_session.vm_glcore ,  
                                                                               name          = 'surface'                 ,
                                                                               active        = True                      ,
                                                                               indexes       = []                        ,
                                                                               is_dynamic    = False                     ,
                                                                               iso_color     = color_minus               ,
                                                                               surface_name  = 'obital_minus'           )
            
            vobject_tmp.parameters = parameters
            
            vobject_tmp.frames = vismol_object.frames
            vobject_tmp.active = False
            vobject_tmp.is_surface = True
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            
            self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
            # Add the VisMol object to the vobject liststore dictionary
            self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
            # Refresh the widgets in the main window
            self.main.refresh_widgets()
            self.vm_session.vm_glcore.queue_draw()
            self.counter +=1

        else:
            pass
        
        
        if vobject_tmp:
            #setting vobject as active in main treeview.
            vobject_tmp.active = True
            self.main.main_treeview.treestore.set_value(vobject_tmp.e_treeview_iter, 6, True)





    def on_button_import_wavefunction (self, widget):
        """ Function doc """
        print('on_button_import_wavefunction')
        
        system_id = self.system_names_combo.get_system_id()
        system    = self.main.p_session.psystem[system_id]
        
        vobject_id    = self.coordinates_combobox.get_vobject_id()
        vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]
        print(system_id, vobject_id, system, vismol_object)
        backup = []
        try:
            backup.append(system.e_treeview_iter)
            backup.append(system.e_liststore_iter)
            system.e_treeview_iter   = None
            system.e_liststore_iter  = None
        except:
            pass
        
        
        #frame = self.vm_session.frame
        #coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
        #joblist = [[frame, system, coords]] 
        
        #'''
        trajectory = [None]*vismol_object.frames.shape[0]
        joblist = []
        for frame in range(vismol_object.frames.shape[0]):
            self.p_session.set_psystem_coordinates_from_vobject( vobject = vismol_object, 
                                                                          system_id = system_id, 
                                                                          frame = frame)
            parameters = None
            coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
            joblist.append([frame, system, coords])
        #'''    
        p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
        results = p.map(generate_wavefunction_parallel, joblist)
        
        self.wave_function_dict[vobject_id] = results
        

        try:
            system.e_treeview_iter   = backup[0]
            system.e_liststore_iter  = backup[1]
        except:
            pass

        #print(self.wave_function_dict)


        
        #'''
        
        self.orbital_liststore_dict[vobject_id]= []
        
        for frame , data in enumerate(self.wave_function_dict[vobject_id]):
            orbitals = data[0]
            
            self.liststore = Gtk.ListStore(int, str, int, float, bool)
            for i in range(len(orbitals)):
                reverse_index = -i-1 #- len(orbitals)
                #print(reverse_index, orbitals[reverse_index ])
        
                self.liststore.append(orbitals[reverse_index ])
                self.orbital_liststore_dict[vobject_id].append(self.liststore)
        
        print()
        self.treeview.set_model(self.orbital_liststore_dict[vobject_id][self.frame])


    def on_treeview_Objects_button_release_event(self, tree, event):
        '''
         str  ,   #                                   # 0
         bool ,   # toggle active=1                   # 1
         bool ,   # toggle visible = 3                # 2 
                                                      
         bool ,   # radio  active  = 2                # 3 
         bool ,   # radio  visible = 4                # 4 
                                                      
         bool  ,  # traj radio  active = 5            # 5 
         bool  ,  # is trajectory radio visible?      # 6 
                                                      
         int,     #                                   # 7
         int,     # pdynamo system index              # 8
         int,)    # frames  # 9
        '''
        
        
        _id = self.system_names_combo.get_active()
        if _id == -1:
            '''_id = -1 means no item inside the combobox'''
            return None
        else:    
            _, system_id = self.main.system_liststore[_id]
            
            
            
        if event.button == 3:
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
            for item in model:
                pass
                #print (item[0], model[iter][0])
            if iter != None:
                self.treeview_menu.open_menu(iter, system_id)

        if event.button == 1:
            print ('event.button == 1')



# NOTA: compute_smooth_normals() (calculo de normal por vertice via media
# das normais de face, em Python puro) foi removida -- substituida pelo
# metodo nativo surface.MakeVertexNormalsFromPolygonalNormals(), que ja
# vem com o pDynamo3 (PolygonalSurface.py) e implementa exatamente o
# mesmo algoritmo, so que compilado (Cython/C) em vez de loop Python.
# Ver surface_parser() e surface_parser_mep() logo abaixo.


def cube_to_pdynamo_surface ( cube_grid, isovalue ):
    # [EN] Bridges an externally-read .cube file (util.cube_reader.CubeGrid)
    # into pDynamo3's own compiled marching-cubes routine, so that a file
    # from ORCA/orca_plot (or any other Gaussian-Cube-compatible program)
    # can be turned into a PolygonalSurface using the SAME code path pDynamo3
    # already uses for its own QC-system-derived orbitals/density/potential
    # surfaces. See item 9 of the changelog at the top of this file for the
    # full story (why RegularGrid.FromDimensionData + Array.FromIterable +
    # Reshape + MarchingCubes_Isosurface3D was chosen over e.g. scikit-image,
    # and why this could not be executed in the assistant's own environment).
    """ Recebe um CubeGrid (de util.cube_reader.read_cube_file) e um
    isovalor, devolve um objeto PolygonalSurface NATIVO do pDynamo3 --
    o mesmo tipo de objeto que generator.GetProperty(tag).isosurface ja
    devolve pros outros tipos de superficie (orbital/potential/density),
    pronto pra passar direto em surface_parser()/surface_parser_mep()
    sem nenhuma adaptacao adicional.

    Reusa o marching cubes compilado do proprio pDynamo3
    (MarchingCubes_Isosurface3D), que so pede um RegularGrid + um array
    de valores + um isovalor -- nenhuma dependencia de sistema QC do
    pDynamo. Ver pScientific/Surfaces/__extensions__/pyrex/
    pScientific.Surfaces.MarchingCubes.pyx e pSimulation/QCGridProperties.py
    no codigo-fonte do pDynamo3 (github.com/pdynamo/pDynamo3). """
    if not cube_grid.is_orthogonal:
        raise ValueError (
            "RegularGrid do pDynamo so aceita grids com eixos alinhados "
            "(vetores de voxel diagonais, sem rotacao/cisalhamento). Este "
            ".cube tem eixos nao-ortogonais -- incomum, mas alguns programas "
            "podem gerar isso; nao suportado por este importador."
        )
    nx, ny, nz = cube_grid.dims
    dx, dy, dz = cube_grid.spacing
    ox, oy, oz = cube_grid.origin

    grid = RegularGrid.FromDimensionData ( [
        { "bins" : nx, "binSize" : dx, "lower" : ox } ,
        { "bins" : ny, "binSize" : dy, "lower" : oy } ,
        { "bins" : nz, "binSize" : dz, "lower" : oz } ,
    ] )

    # cube_grid.values ja esta na ordem (nx,ny,nz) com Z mais rapido
    # (ver cube_reader.py) -- mesma convencao de reshape que
    # QCGridProperty.Isosurface() usa internamente pro grid do proprio
    # pDynamo, entao um reshape simples em ordem C basta.
    flat_values = cube_grid.values.reshape ( -1 ).tolist ( )
    flat_array  = Array.FromIterable ( flat_values )
    dataND      = Reshape ( flat_array, ( nx, ny, nz ), resultClass = RealArrayND )

    surface = MarchingCubes_Isosurface3D ( grid, dataND, isovalue )
    surface.MakeVertexNormalsFromPolygonalNormals ( )
    return surface


def _trilinear_interpolate ( values_3d, origin, spacing, query_points ):
    """ [EN] Hand-written, pure-numpy trilinear interpolation over a regular
    3D grid -- no scipy needed. Replaces the previous
    scipy.interpolate.RegularGridInterpolator (for .cube grids) and
    scipy.spatial.cKDTree nearest-neighbour lookup (for pDynamo grids,
    which used to be necessary only because we didn't yet know pDynamo's
    RegularGrid exposes .origin/.spacing/.shape directly -- see
    build_potential_interpolator() below).

    values_3d : array (nx, ny, nz)      -- scalar field on the grid
    origin    : (ox, oy, oz)            -- world coordinate of grid index (0,0,0)
    spacing   : (dx, dy, dz)            -- grid spacing along each axis
    query_points : array (m, 3)         -- world-space points to evaluate,
                                            same units as origin/spacing (Bohr)

    Devolve um array (m,) com os valores interpolados. Pontos fora da caixa
    do grid sao GRAMPEADOS (clipped) pro voxel mais proximo da borda, em
    vez de extrapolar ou quebrar -- comportamento razoavel aqui, ja que a
    caixa do grid normalmente ja envolve toda a superficie de interesse. """
    nx, ny, nz = values_3d.shape
    ox, oy, oz = origin
    dx, dy, dz = spacing
    query_points = np.asarray ( query_points, dtype = np.float64 )

    # world coords -> fractional grid-index coords
    fx = ( query_points[:,0] - ox ) / dx
    fy = ( query_points[:,1] - oy ) / dy
    fz = ( query_points[:,2] - oz ) / dz

    # clamp so that i0+1 never goes out of bounds (handles points slightly
    # outside the box, e.g. a density-isosurface vertex sitting right at
    # the edge of the potential grid's own bounding box)
    fx = np.clip ( fx, 0.0, nx - 1 - 1e-9 )
    fy = np.clip ( fy, 0.0, ny - 1 - 1e-9 )
    fz = np.clip ( fz, 0.0, nz - 1 - 1e-9 )

    ix0 = np.floor ( fx ).astype ( np.int64 ); ix1 = ix0 + 1
    iy0 = np.floor ( fy ).astype ( np.int64 ); iy1 = iy0 + 1
    iz0 = np.floor ( fz ).astype ( np.int64 ); iz1 = iz0 + 1

    tx = fx - ix0
    ty = fy - iy0
    tz = fz - iz0

    # the 8 corners of the voxel containing each query point
    c000 = values_3d[ix0, iy0, iz0]; c001 = values_3d[ix0, iy0, iz1]
    c010 = values_3d[ix0, iy1, iz0]; c011 = values_3d[ix0, iy1, iz1]
    c100 = values_3d[ix1, iy0, iz0]; c101 = values_3d[ix1, iy0, iz1]
    c110 = values_3d[ix1, iy1, iz0]; c111 = values_3d[ix1, iy1, iz1]

    # interpolate along x, then y, then z (standard trilinear recipe)
    c00 = c000 * ( 1 - tx ) + c100 * tx
    c01 = c001 * ( 1 - tx ) + c101 * tx
    c10 = c010 * ( 1 - tx ) + c110 * tx
    c11 = c011 * ( 1 - tx ) + c111 * tx

    c0 = c00 * ( 1 - ty ) + c10 * ty
    c1 = c01 * ( 1 - ty ) + c11 * ty

    return c0 * ( 1 - tz ) + c1 * tz


def build_potential_interpolator_from_cube ( cube_grid ):
    # [EN] External-.cube analogue of build_potential_interpolator() below.
    # A parsed CubeGrid already gives an explicit regular grid (known
    # origin + per-axis spacing + values already shaped (nx,ny,nz)), so
    # this is a thin wrapper around _trilinear_interpolate(). See
    # changelog item 9 for the cube-import feature this supports.
    """ Como build_potential_interpolator(), mas para um CubeGrid lido de
    um arquivo .cube externo (ORCA, etc) em vez de um QCGridProperty do
    pDynamo. """
    nx, ny, nz = cube_grid.dims
    dx, dy, dz = cube_grid.spacing
    ox, oy, oz = cube_grid.origin

    def evaluate ( query_points ):
        """ query_points: array (m,3), nas MESMAS unidades do grid (Bohr). """
        return _trilinear_interpolate ( cube_grid.values, (ox,oy,oz), (dx,dy,dz), query_points )

    return evaluate


def _is_degenerate_vertex ( vertices, v ):
    """ Deteccao rapida e barata (mas incompleta sozinha -- ver
    _compute_valid_polygon_mask) de um vertice "fantasma" que caiu
    exatamente na origem do MUNDO (0.0, 0.0, 0.0). So cobre o caso em
    que a origem do grid (ver nota abaixo) coincide com a origem do
    mundo -- o que acontece as vezes (ex: nosso .cube sintetico de
    teste, centrado propositalmente na origem), mas nao em geral. """
    return ( vertices[v,0] == 0.0 ) and ( vertices[v,1] == 0.0 ) and ( vertices[v,2] == 0.0 )


def _compute_valid_polygon_mask ( polygons, vertices, size_factor = 8.0 ):
    # [EN] Filters out the marching-cubes "ghost vertex" artifact reported
    # by the user (screenshot showed long spurious triangles shooting off
    # from the surface to one distant point). Root cause (see changelog
    # item 11): unfilled vertex slots in pDynamo3's C marching-cubes buffer
    # keep their zero-initialised raw (grid-index-space) value and, after
    # the library's own final scale+translate step, end up sitting exactly
    # at the REGULAR GRID's own origin corner (RegularGrid's
    # "midPointLower") -- not at the world coordinate origin, which is why
    # an earlier, position-based filter (_is_degenerate_vertex(), kept
    # below purely for reference / no longer called) was insufficient in
    # general. This filter is purely geometric instead: a real marching-
    # cubes triangle is never much bigger than one grid cell, so any
    # triangle whose longest edge is an outlier relative to the mesh's own
    # median edge length is discarded as a ghost-vertex artifact.
    """ Filtro geometrico pros triangulos "fantasma" do marching cubes
    do pDynamo3 (artefato relatado pelo usuario, ver imagem -- linhas
    compridas saindo da superficie ate um ponto isolado longe dela).

    Descoberta lendo MarchingCubes.c do pDynamo3: os vertices sao
    calculados em unidades de INDICE DE GRID (0, 1, 2, ...) e SO NO
    FINAL sao escalados (* binSize) e TRANSLADADOS (+ midPointLower,
    a origem do proprio grid) pra virar coordenada real. Ou seja, um
    vertice "fantasma" que nunca foi tocado pelo algoritmo (fica com o
    valor bruto de inicializacao do array, 0.0 em tudo) NAO acaba na
    origem do MUNDO -- acaba na origem do PROPRIO GRID (midPointLower),
    que pode estar em qualquer canto da caixa, longe da molecula. Por
    isso um filtro checando "vertice == (0,0,0) do mundo" (
    _is_degenerate_vertex acima) so pega o caso em que a origem do grid
    calha de coincidir com a origem do mundo -- nao pega o caso geral.

    A correcao que funciona em geral e geometrica, nao posicional: um
    triangulo de marching cubes de verdade nunca e muito maior que uma
    celula do grid (a malha inteira e local, feita de triangulos
    pequenos e regulares). Um triangulo ligando um vertice bom a um
    vertice fantasma longe tem uma aresta ordens de grandeza maior que
    o normal. Calculamos a maior aresta de cada triangulo, comparamos
    com a MEDIANA de todos os triangulos da malha (auto-calibrado, nao
    precisa saber o espacamento do grid) e descartamos qualquer
    triangulo cuja maior aresta seja muito maior (> size_factor vezes)
    que essa mediana.

    Devolve (mask, n_descartados) -- mask e um array booleano, uma
    entrada por linha de `polygons` (True = triangulo valido, mantem)."""
    n_tri = polygons.rows
    if n_tri == 0:
        return np.zeros ( 0, dtype = bool ), 0

    max_edge = np.empty ( n_tri, dtype = np.float64 )
    for p in range ( n_tri ):
        i0, i1, i2 = polygons[p,0], polygons[p,1], polygons[p,2]
        p0 = np.array ( [ vertices[i0,0], vertices[i0,1], vertices[i0,2] ] )
        p1 = np.array ( [ vertices[i1,0], vertices[i1,1], vertices[i1,2] ] )
        p2 = np.array ( [ vertices[i2,0], vertices[i2,1], vertices[i2,2] ] )
        e01 = np.linalg.norm ( p1 - p0 )
        e12 = np.linalg.norm ( p2 - p1 )
        e20 = np.linalg.norm ( p0 - p2 )
        max_edge[p] = max ( e01, e12, e20 )

    typical = np.median ( max_edge )
    if typical == 0.0:
        typical = 1e-9   # malha degenerada por completo -- evita divisao por zero adiante, sem quebrar
    threshold = size_factor * typical
    mask = max_edge <= threshold
    n_discarded = int ( n_tri - mask.sum ( ) )
    return mask, n_discarded


def surface_parser ( surface, iso_color):
    """ Function doc """
    normals   = surface.polygonNormals
    polygons  = surface.polygons
    vertices  = surface.vertices
    # normal por vertice ja nativa do pDynamo (MarchingCubes.pyx aloca
    # surface.vertexNormals, mas QCGridPropertyGenerator.Isosurface() so
    # chama MakePolygonNormals() -- vertexNormals fica zerado ate a gente
    # chamar isso aqui. Mesmo algoritmo que nosso compute_smooth_normals
    # fazia em Python puro (media das normais de face por vertice,
    # depois normaliza), so que compilado (ver PolygonalSurface.py do
    # pDynamo3).
    surface.MakeVertexNormalsFromPolygonalNormals ( )
    smooth_normals = surface.vertexNormals
    valid_mask, skipped = _compute_valid_polygon_mask ( polygons, vertices )
    colors    = []
    vertices2 = []
    normals2  = []
    
    for p in range ( polygons.rows ):
        #print ( "facet normal " )
        #print ( "\n    outer loop" )
        if not valid_mask[p]:
            continue   # triangulo "fantasma" (aresta anormalmente grande) -- ver _compute_valid_polygon_mask
        tri = polygons[p,:]
        for v in tri:
            #text = "\n        vertex "
            for c in range ( 3 ): 
                vertices2.append((vertices[v,c])/1.889725989 ) # convert from Bohr to angstrom
                #text += " {} ".format(vertices[v,c])
            for c in range ( 3 ):
                normals2.append ( smooth_normals[v, c] )
            for rgb in iso_color: 
                #vertices2.append(rgb)
                colors.append(rgb)
    if skipped:
        print ( "surface_parser: {} triangulo(s) fantasma (aresta anormal) descartado(s)".format ( skipped ) )
    
    vertices = np.array(vertices2, dtype=np.float32)
    colors   = np.array(colors, dtype=np.float32)
    normals  = np.array(normals2, dtype=np.float32)
    # um indice por vertice (nao por componente/float) -- a versao antiga
    # gerava indexes com 3x mais entradas do que vertices de verdade
    # existem no buffer (passava despercebido, mas era um out-of-bounds
    # read em potencial na GPU -- ver nota no README).
    indexes  = np.array(range(len(vertices)//3), dtype=np.uint32)
    return vertices, colors, indexes, normals


def _colormap_lookup ( t_values, color_map ):
    """ Interpola um colormap no formato de util/colormaps.py (dict
    posicao [0..~1] -> [r,g,b]) de forma vetorizada com numpy.interp --
    equivalente ao get_color()/interpolate_color() de
    util/easyplot/color_utils.py, mas para um array inteiro de uma vez
    (mais rapido pra milhares de vertices, e evita importar color_utils.py
    -- que arrasta gi/Gtk/cairo so pra reusar uma funcao escalar -- dentro
    dos processos worker do multiprocessing.Pool). """
    stops  = np.array ( sorted ( color_map.keys() ), dtype = np.float64 )
    colors = np.array ( [ color_map[s] for s in stops ], dtype = np.float64 )
    r = np.interp ( t_values, stops, colors[:, 0] )
    g = np.interp ( t_values, stops, colors[:, 1] )
    b = np.interp ( t_values, stops, colors[:, 2] )
    return np.stack ( [r, g, b], axis = 1 ).astype ( np.float32 )


def mep_colormap ( values, vmin = None, vmax = None, cmap_name = 'coolwarm', reverse = True, percentile = 2.0 ):
    # [EN] Maps a 1-D array of per-vertex electrostatic-potential values to
    # RGB colour, for the MEP (Molecular Electrostatic Potential) surface
    # type -- see changelog item 6 for the full history (matplotlib was
    # tried first, then dropped in favour of the project's own COLOR_MAPS
    # dict, see _colormap_lookup() just above). The `percentile` parameter
    # (added after live debugging, see changelog item 8) is important: a
    # single outlier vertex (e.g. one sitting extremely close to a nucleus,
    # where the potential genuinely diverges) can otherwise single-handedly
    # define the whole colour scale via naive min()/max(), squashing every
    # other, "normal" vertex to within a fraction of a percent of the exact
    # middle of the range -- i.e. a uniform, washed-out grey instead of a
    # gradient.
    """ Mapeia valores escalares (potencial eletrostatico por vertice) para
    RGB usando os colormaps definidos em util/colormaps.py (COLOR_MAPS),
    em vez dos colormaps do matplotlib.

    cmap_name: qualquer chave de COLOR_MAPS -- 'coolwarm' (default), 'vik',
    'berlin', 'bam  ' (com os 2 espacos, e o nome literal no dict), 'jet',
    'rainbow', 'gnuplot', 'magma', 'viridis', 'plasma', 'cividis', etc.
    Ver util/colormaps.py pra lista completa e pontos de controle de cada um.

    O centro (potencial = 0) e ancorado em t = 0.5 mesmo quando vmin/vmax
    nao sao simetricos -- equivalente ao TwoSlopeNorm que usavamos com
    matplotlib.colors, so que implementado a mao aqui.

    percentile (default 2.0): se vmin/vmax NAO forem passados
    explicitamente, em vez de usar o minimo/maximo brutos dos valores,
    usa os percentis [percentile, 100-percentile]. Isso importa porque
    um unico vertice patologico -- por exemplo, muito perto de um nucleo,
    onde o potencial eletrostatico diverge como 1/r -- pode ter uma
    magnitude ordens de grandeza maior que o resto da superficie (ja visto
    na pratica: -0.25 a.u. tipico vs. 54000+ a.u. num vertice isolado). Se
    isso vira o limite da escala, TODO o resto da superficie fica
    esmagado a menos de 0.01% de distancia do centro -- aparece como uma
    cor quase uniforme (cinza claro / "branco"). Com o corte por
    percentil, esses poucos vertices extremos so ficam saturados na cor
    mais forte de uma ponta (clip do np.clip abaixo), sem distorcer a
    escala de todo mundo.

    reverse=True (default): os mapas divergentes de COLOR_MAPS (coolwarm,
    vik, berlin) vao de azul (posicao 0.0) pra vermelho (posicao ~1.0) --
    o oposto da convencao quimica que queremos pro MEP (vermelho =
    negativo/regiao rica em eletrons, azul = positivo/regiao pobre em
    eletrons). reverse=True inverte a consulta (t_query = 1 - t) pra
    bater com essa convencao. Pra mapas sequenciais tipo 'jet'/'rainbow'
    (sem centro definido), reverse ainda funciona mas a nocao de "centro
    em zero" fica menos significativa -- nesse caso convem chamar com
    vmin/vmax explicitos cobrindo o range real dos dados. """
    values = np.asarray ( values, dtype = np.float64 )
    vmin = float ( vmin if vmin is not None else np.percentile ( values, percentile ) )
    vmax = float ( vmax if vmax is not None else np.percentile ( values, 100.0 - percentile ) )

    if vmin >= 0:
        vmin = -abs ( vmax ) if vmax != 0 else -1e-9
    if vmax <= 0:
        vmax = abs ( vmin ) if vmin != 0 else 1e-9
    limit = max ( abs ( vmin ), abs ( vmax ) )
    if limit == 0:
        limit = 1e-9

    # t = 0.5 no zero, 0.0 no -limit, 1.0 no +limit -- simetrico em torno
    # do centro (o lado "mais estreito" da faixa so nunca chega a atingir
    # t = 0 ou t = 1 exatos, o que e o comportamento correto). Valores
    # além de +-limit (os outliers cortados pelo percentil) sao
    # grampeados em t=0 ou t=1 pelo np.clip -- ficam saturados, nao
    # quebram a escala.
    t = 0.5 + 0.5 * np.clip ( values / limit, -1.0, 1.0 )

    if reverse:
        t = 1.0 - t

    return _colormap_lookup ( t, COLOR_MAPS[cmap_name] )


def _nearest_neighbor_lookup ( pts, vals, query_points, chunk_size = 4000 ):
    """ [EN] Pure-numpy replacement for scipy.spatial.cKDTree, used ONLY for
    pDynamo3's own potential grid (build_potential_interpolator() below) --
    NOT for external .cube files, which use real trilinear interpolation
    instead (_trilinear_interpolate() above), because for a .cube file we
    control the parsing ourselves and know its (nx,ny,nz) axis order for
    certain (see util/cube_reader.py).

    pDynamo3's grid, by contrast, only exposes gridPoints/gridValues as two
    flat PARALLEL arrays (point <-> value by matching row index) -- we do
    NOT actually know, confirmed, whether reshaping gridValues straight
    into potentialProperty.grid.shape and treating it as a simple C-order
    (x slowest, z fastest) 3D array would land on the same convention
    pDynamo3 uses internally for its own RegularGrid. Guessing wrong there
    would silently produce a WRONG but plausible-looking interpolated
    surface -- not a crash, just quietly incorrect chemistry. Nearest-
    neighbour sidesteps the whole question: it works directly off the
    flat (point, value) pairs, with no assumption about their storage
    order at all, at the cost of being a little less accurate than true
    trilinear (acceptable here, given typical grid spacings).

    pts  : array (n, 3) -- grid point coordinates, any order
    vals : array (n,)   -- matching scalar values, same order as pts
    query_points : array (m, 3)
    chunk_size: caps memory use of the (chunk, n, 3) distance tensor built
    per batch -- with n up to ~150000 (a fairly fine QC grid) and
    chunk_size=4000, that is at most 4000*150000*3 floats (~7 GB worst
    case if done in one shot without chunking; chunking keeps peak memory
    to chunk_size*n*3 floats at a time, a few hundred MB, and runs in a
    handful of vectorised numpy passes rather than one huge allocation). """
    pts = np.asarray ( pts, dtype = np.float64 )
    vals = np.asarray ( vals, dtype = np.float64 )
    query_points = np.asarray ( query_points, dtype = np.float64 )
    m = query_points.shape[0]
    out = np.empty ( m, dtype = np.float64 )
    for start in range ( 0, m, chunk_size ):
        end = min ( start + chunk_size, m )
        chunk = query_points[start:end]                       # (c,3)
        diff = chunk[:, None, :] - pts[None, :, :]             # (c,n,3)
        d2 = np.einsum ( 'ijk,ijk->ij', diff, diff )           # (c,n) squared distances
        idx = np.argmin ( d2, axis = 1 )                       # (c,) index of nearest grid point
        out[start:end] = vals[idx]
    return out


def build_potential_interpolator ( potentialProperty ):
    # [EN] Uses nearest-neighbour (_nearest_neighbor_lookup() above), NOT
    # trilinear interpolation, and NOT scipy.spatial.cKDTree either --
    # see the long comment on _nearest_neighbor_lookup() for exactly why:
    # in short, pDynamo3's own grid object only gives flat, parallel
    # (point, value) arrays with no confirmed/known [i,j,k] storage order,
    # so we cannot safely reshape gridValues into (nx,ny,nz) and reuse the
    # exact-and-verified _trilinear_interpolate() path the way
    # build_potential_interpolator_from_cube() does for external .cube
    # files (where WE control the parsing and the order is known for
    # certain). This removes the scipy dependency without introducing an
    # unverified assumption about pDynamo3's internal array layout.
    """ Recebe o QCGridProperty bruto do potencial (generator.GetProperty(tag),
    ANTES de virar isosuperficie -- precisa ter .gridPoints/.gridValues) e
    devolve uma funcao que avalia o potencial em qualquer ponto 3D via
    vizinho mais proximo no grid denso do pDynamo. """
    n = len ( potentialProperty.gridValues )
    pts = np.empty ( (n, 3), dtype = np.float64 )
    for i in range ( n ):
        pts[i, 0] = potentialProperty.gridPoints[i, 0]
        pts[i, 1] = potentialProperty.gridPoints[i, 1]
        pts[i, 2] = potentialProperty.gridPoints[i, 2]
    vals = np.array ( [ potentialProperty.gridValues[i] for i in range ( n ) ], dtype = np.float64 )

    def evaluate ( query_points ):
        """ query_points: array (m,3), nas MESMAS unidades do grid (Bohr,
        que e a unidade nativa do pDynamo -- ver surface_parser, que so
        converte pra Angstrom na hora de montar o buffer de exibicao). """
        return _nearest_neighbor_lookup ( pts, vals, query_points )

    return evaluate


def surface_parser_mep ( surface, vertex_colors ):
    """ Como surface_parser, mas recebe uma cor RGB ja calculada por vertice
    (vertex_colors, shape (n_vertices_originais, 3)) em vez de um iso_color
    unico repetido pra malha inteira -- usado pro MEP (mapa continuo de
    potencial eletrostatico sobre a superficie de densidade). """
    polygons  = surface.polygons
    vertices  = surface.vertices
    # ver nota em surface_parser() -- usa o metodo nativo do pDynamo em
    # vez do compute_smooth_normals em Python puro.
    surface.MakeVertexNormalsFromPolygonalNormals ( )
    smooth_normals = surface.vertexNormals
    valid_mask, skipped = _compute_valid_polygon_mask ( polygons, vertices )
    colors    = []
    vertices2 = []
    normals2  = []

    for p in range ( polygons.rows ):
        if not valid_mask[p]:
            continue   # triangulo "fantasma" (aresta anormalmente grande) -- ver _compute_valid_polygon_mask
        tri = polygons[p, :]
        for v in tri:
            for c in range ( 3 ):
                vertices2.append ( (vertices[v, c]) / 1.889725989 )  # Bohr -> Angstrom
            for c in range ( 3 ):
                normals2.append ( smooth_normals[v, c] )
            for rgb in vertex_colors[v]:
                colors.append ( rgb )
    if skipped:
        print ( "surface_parser_mep: {} triangulo(s) fantasma (aresta anormal) descartado(s)".format ( skipped ) )

    vertices = np.array ( vertices2, dtype = np.float32 )
    colors   = np.array ( colors, dtype = np.float32 )
    normals  = np.array ( normals2, dtype = np.float32 )
    # um indice por vertice (nao por componente/float) -- ao contrario do
    # surface_parser original, que gera indexes com 3x mais entradas do
    # que vertices de verdade existem no buffer (ver nota no README).
    indexes  = np.array ( range ( len(vertices) // 3 ), dtype = np.uint32 )
    return vertices, colors, indexes, normals


def apply_coords_to_system (system, coords):
    """ Function doc """
    
    for i, xyz in enumerate(coords):
        system.coordinates3[i][0] = xyz[0]
        system.coordinates3[i][1] = xyz[1]
        system.coordinates3[i][2] = xyz[2]


def _generate_external_cube_surface ( parameters ):
    # [EN] Entry point for the "External" surface type (changelog item 9):
    # generates a surface entirely from user-supplied .cube file(s) (e.g.
    # exported from ORCA via orca_plot), with NO live pDynamo QC system
    # involved at all -- called from an early-return guard at the very top
    # of generate_grid_parallel(), before that function's normal,
    # unconditional QC-system setup (apply_coords_to_system / system.Energy
    # / QCGridPropertyGenerator.FromSystem) would otherwise run and crash on
    # a None system/coords pair.
    """ Gera a superficie a partir de um arquivo .cube externo (ex: ORCA
    via orca_plot), sem nenhuma dependencia de sistema QC do pDynamo --
    so leitura de arquivo + marching cubes nativo (ver
    cube_to_pdynamo_surface / util.cube_reader). """
    density_path   = parameters['external_density_path']
    potential_path = parameters.get ( 'external_potential_path' )
    isovalue       = parameters['_isovalue']
    color_plus     = parameters['color_plus']
    mep_vmin       = parameters.get ( 'mep_vmin' )
    mep_vmax       = parameters.get ( 'mep_vmax' )
    mep_cmap_name  = parameters.get ( 'mep_cmap_name', 'coolwarm' )

    density_cube = read_cube_file ( density_path )
    surface      = cube_to_pdynamo_surface ( density_cube, isovalue )

    if potential_path:
        # colore por potencial eletrostatico interpolado (MEP), igual ao
        # branch 'mep' -- so que os dois cubos (densidade e potencial)
        # vem de arquivos externos em vez do grid do pDynamo.
        potential_cube     = read_cube_file ( potential_path )
        evaluate_potential = build_potential_interpolator_from_cube ( potential_cube )
        n_verts    = surface.vertices.rows
        verts_bohr = np.empty ( (n_verts, 3), dtype = np.float64 )
        for v in range ( n_verts ):
            for c in range ( 3 ):
                verts_bohr[v, c] = surface.vertices[v, c]
        potential_values = evaluate_potential ( verts_bohr )
        vertex_colors    = mep_colormap ( potential_values, vmin = mep_vmin, vmax = mep_vmax,
                                           cmap_name = mep_cmap_name )
        vertices, colors, indexes, normals = surface_parser_mep ( surface, vertex_colors )
    else:
        # sem cubo de potencial -- cor uniforme, igual density/orbital.
        vertices, colors, indexes, normals = surface_parser ( surface, iso_color = color_plus )

    return { 'obital_plus' : [ vertices, colors, indexes, normals ] }


def generate_grid_parallel (job):
    """ Function doc 
    
    [frame, system, coords, parameters]
    
    """
    i          = job[0]
    system     = job[1]
    coords     = job[2]
    parameters = job[3]
    #_type      = job[4]
    '''
    parameters = {
            '_GridSpacing'   : _GridSpacing
            '_OrbitalTag'    : _OrbitalTag
            '_IsosurfaceTag' : _IsosurfaceTag
            }
    '''
    _type = parameters['type']

    # [EN] Early-return guard, added for the "External" cube-import surface
    # type (changelog item 9). Every other branch below unconditionally
    # calls apply_coords_to_system() / system.Energy() /
    # QCGridPropertyGenerator.FromSystem(system) a few lines further down,
    # all of which require a real, live pDynamo QC system -- but an
    # externally-supplied .cube file has no such system (system/coords are
    # simply None in that job tuple). Must return BEFORE that setup code
    # runs, not after.
    if _type == 'external_cube':
        # cubo externo (.cube, ex: ORCA via orca_plot) -- nao precisa de
        # system/coords/QCGridPropertyGenerator nenhum, so leitura de
        # arquivo. Desvia ANTES do setup de sistema QC logo abaixo, que
        # exigiria um "system" de verdade (None aqui, ja que essa entrada
        # nao vem de um calculo QC do pDynamo).
        return _generate_external_cube_surface ( parameters )

    _GridSpacing   = parameters['_GridSpacing']
    _OrbitalTag    = parameters['_OrbitalTag']
    _isovalue      = parameters['_isovalue']  
    _IsosurfaceTag = parameters['_IsosurfaceTag']
    _mep_vmin      = parameters.get ( 'mep_vmin' )   # None = automatico (percentil)
    _mep_vmax      = parameters.get ( 'mep_vmax' )
    _mep_cmap_name = parameters.get ( 'mep_cmap_name', 'coolwarm' )
    key            = parameters['orbital_key']
    color_plus     = parameters['color_plus']
    color_minus    = parameters['color_minus']
    
    #print(parameters, type(system))
    #-----------------------------------------------------------------------
    # . Calculate the system grid properties.
    #-----------------------------------------------------------------------
    #system    = system.Energy()
    #energies  = orbitalsP.energies
    
    apply_coords_to_system(system, coords)
    system.Energy()
    
    #-----------------------------------------------------------------------
    # . Calculate the system grid properties.
    #-----------------------------------------------------------------------
    generator = QCGridPropertyGenerator.FromSystem (system )
    generator.DefineGrid    ( gridSpacing = _GridSpacing ) # . Some value in atomic units - e.g. 0.2
    

    
    orbital_iso = {}
    if _type == 'orbital':
        generator.GridOrbitals  ( [ key ]    ,       tag = _OrbitalTag    ) # . List of orbital indices (can be one only)    
        generator.Isosurface    ( _OrbitalTag, _isovalue, tag = _IsosurfaceTag )
    
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_p = surfaceProperty.isosurface # . This is the surface you can display.
        
        vertices, colors, indexes, normals = surface_parser ( surface = isosurface_p , iso_color = color_plus )
        
        orbital_iso['obital_plus'] = [vertices, colors, indexes, normals]
        generator.Isosurface    ( _OrbitalTag, _isovalue*-1, tag = _IsosurfaceTag )
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_n = surfaceProperty.isosurface # . This is the surface you can display.
        
        vertices, colors, indexes, normals = surface_parser ( surface = isosurface_n , iso_color = color_minus )
        orbital_iso['obital_minus'] = [vertices, colors, indexes, normals]
    
    elif _type == 'potential':
        generator.GridPotential  (                   tag = 'potential'      ) # . List of orbital indices (can be one only)    
        generator.Isosurface    ( 'potential', _isovalue, tag = _IsosurfaceTag)
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_p = surfaceProperty.isosurface # . This is the surface you can display.
        
        vertices, colors, indexes, normals = surface_parser ( surface = isosurface_p , iso_color = color_plus )
        
        orbital_iso['obital_plus'] = [vertices, colors, indexes, normals]
        generator.Isosurface    ( 'potential', _isovalue*-1, tag = _IsosurfaceTag )
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_n = surfaceProperty.isosurface # . This is the surface you can display.
        
        vertices, colors, indexes, normals = surface_parser ( surface = isosurface_n , iso_color = color_minus )
        orbital_iso['obital_minus'] = [vertices, colors, indexes, normals]
    
    elif _type == 'mep':
        # 1. Geometria da malha a partir da isosuperficie de DENSIDADE
        #    (o isovalor do campo entry_isovalue passa a significar
        #    "isovalor de densidade" nesse modo -- ~0.002-0.02 u.a. costuma
        #    aproximar bem o contorno de van der Waals).
        generator.GridDensity ( tag = 'density_mep' )
        generator.Isosurface  ( 'density_mep', _isovalue, tag = _IsosurfaceTag )
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        density_iso = surfaceProperty.isosurface

        # 2. Grid de POTENCIAL bruto (os valores a mapear na malha acima).
        #    Tag propria, distinta de 'density_mep' e de _IsosurfaceTag --
        #    ver o bug de colisao de tags que corrigimos no branch 'potential'.
        generator.GridPotential ( tag = 'potential_mep' )
        potentialProperty  = generator.GetProperty ( 'potential_mep' )
        evaluate_potential = build_potential_interpolator ( potentialProperty )

        # 3. Avalia o potencial em cada vertice ORIGINAL da malha de densidade,
        #    em Bohr (mesma unidade do grid do pDynamo -- a conversao pra
        #    Angstrom so acontece dentro de surface_parser_mep).
        n_verts    = density_iso.vertices.rows
        verts_bohr = np.empty ( (n_verts, 3), dtype = np.float64 )
        for v in range ( n_verts ):
            for c in range ( 3 ):
                verts_bohr[v, c] = density_iso.vertices[v, c]

        potential_values = evaluate_potential ( verts_bohr )
        vertex_colors    = mep_colormap ( potential_values, vmin = _mep_vmin, vmax = _mep_vmax,
                                           cmap_name = _mep_cmap_name )

        vertices, colors, indexes, normals = surface_parser_mep ( density_iso, vertex_colors )
        orbital_iso['obital_plus'] = [vertices, colors, indexes, normals]
    
    
    
    
    
    elif _type == 'density':
        generator.GridDensity  (                   tag = 'density'   ) # . List of orbital indices (can be one only)    
        generator.Isosurface    ( 'density', _isovalue, tag = _IsosurfaceTag)
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_p = surfaceProperty.isosurface # . This is the surface you can display.
        vertices, colors, indexes, normals = surface_parser ( surface = isosurface_p , iso_color = color_plus )
        orbital_iso['obital_plus'] = [vertices, colors, indexes, normals]
        
    else:
        pass
        
   
    #vertices, colors, indexes = surface_parser ( surface = isosurface_p , iso_color = [1,0,0] )
    #
    #orbital_iso['obital_plus'] = [vertices, colors, indexes]
    #
    #generator.Isosurface    ( _OrbitalTag, _isovalue*-1, tag = _IsosurfaceTag )
    #surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
    #isosurface_n = surfaceProperty.isosurface # . This is the surface you can display.
    #
    #vertices, colors, indexes = surface_parser ( surface = isosurface_n , iso_color = [0,0,1] )
    #orbital_iso['obital_minus'] = [vertices, colors, indexes]
    
    return orbital_iso
    
    #generator.DefineGrid    ( gridSpacing = _GridSpacing ) # . Some value in atomic units - e.g. 0.2
    #
    #orbital_iso = {}
    #print ('key, _OrbitalTag:', key, _OrbitalTag)
    #generator.GridOrbitals  ( [ key ], tag = _OrbitalTag) # . List of orbital indices (can be one only)    
    #
    #generator.Isosurface    ( _OrbitalTag, _isovalue, tag = _IsosurfaceTag )
    #surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
    #isosurface_p = surfaceProperty.isosurface # . This is the surface you can display.
    #
    #
    #
    #vertices, colors, indexes = surface_parser ( surface = isosurface_p , iso_color = [1,0,0] )
    #
    #orbital_iso['obital_plus'] = [vertices, colors, indexes]
    #
    #generator.Isosurface    ( _OrbitalTag, _isovalue*-1, tag = _IsosurfaceTag )
    #surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
    #isosurface_n = surfaceProperty.isosurface # . This is the surface you can display.
    #
    #vertices, colors, indexes = surface_parser ( surface = isosurface_n , iso_color = [0,0,1] )
    #orbital_iso['obital_minus'] = [vertices, colors, indexes]
    #
    #return orbital_iso



def generate_wavefunction_parallel(job):
    """ Function doc """
    generate_grid_parallel

    i          = job[0]
    system     = job[1]
    coords     = job[2]
    #parameters = job[3]
    
    apply_coords_to_system(system, coords)
    system.Energy()

    orbitalsP = system.scratch.orbitalsP
    energies  = orbitalsP.energies
    
    LUMO      = orbitalsP.occupancyHandler.numberOccupied
    HOMO      = LUMO - 1
    
    #generator = QCGridPropertyGenerator.FromSystem (system )
    #generator.orbitalsP = orbitalsP
    orbitals  = []
    generator = QCGridPropertyGenerator.FromSystem ( system )
    
    for i,energy in enumerate(energies):
        if i >= LUMO:
            if i-LUMO == 0:
                label = 'LUMO ' 
            else:
                label = 'LUMO +'+str(i-LUMO)
        else:
            
            if i-HOMO == 0:
                label = 'HOMO' 
            else:
                label = 'HOMO '+str(i-HOMO)
        orbitals.append([i, label, orbitalsP.occupancies[i], energy, False ])
        
    return orbitals, system, generator
    
    
    
    
    
    '''
    for i in range(len(orbitals)):
        reverse_index = -i-1 #- len(orbitals)
        print(reverse_index, orbitals[reverse_index ])
        self.liststore.append(orbitals[reverse_index ])
    '''











class TreeViewMenu:
    """ Class doc """
    
    def __init__ (self, sele_window):
        """ Class initialiser """
        pass
        self.treeview = sele_window.treeview
        self.p_session = sele_window.p_session
        self.sele_window = sele_window 
        functions = {
                    'Rename'                : self.print_test ,
                    'Delete'                : self.delete_system ,
                    }
        self.build_tree_view_menu(functions)
        self.rename_window_visible = False

    def print_test (self, menu_item = None ):
        """  
        menu_item = Gtk.MenuItem object at 0x7fbdcc035700 (GtkMenuItem at 0x37cf6c0)
        
        """
        if self.rename_window_visible:
            pass
        else:
            #
            self.e_id     = self.sele_window.system_names_combo.get_active()
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
            self.key      = model.get_value(iter, 0)
            sys           = model.get_value(iter, 1)
            
            #print('key: ', self.key, 'e_id: ',self.e_id)
            self.window = Gtk.Window()
            self.window.connect('destroy', self.destroy)
            self.window.set_keep_above(True)
            self.entry  = Gtk.Entry()
            
            self.entry.connect('activate', self.rename)
            self.window.add(self.entry)
            self.rename_window_visible = True
            self.window.show_all()
            #print(menu_item)

    def rename (self, menu_item):
        """ Function doc """
        #print('New name: ', self.entry.get_text())
        new_name = self.entry.get_text()
        pass
        
        
        self.p_session.psystem[self.e_id].e_selections[new_name] = self.p_session.psystem[self.e_id].e_selections[self.key]
        self.p_session.psystem[self.e_id].e_selections.pop(self.key)
        self.sele_window.update_window()
        self.window.destroy()
        self.rename_window_visible = False
        
    def destroy (self, widget):
        """ Function doc """
        self.rename_window_visible = False

    def delete_system (self, menu_item = None ):
        """ Function doc """
        selection = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        #print(model[iter][0])


        sele = self.p_session.psystem[self.system_id].e_selections.pop(model[iter][0])
        #print ('deleting',sele)
        #print ('selections', self.p_session.systems[self.system_id]['selections'])
        self.sele_window.update_window (system_names = False, coordinates = False,  selections = True )


    def build_tree_view_menu (self, menu_items = None):
        """ Function doc """
        self.tree_view_menu = Gtk.Menu()
        for label in menu_items:
            mitem = Gtk.MenuItem(label)
            mitem.connect('activate', menu_items[label])
            self.tree_view_menu.append(mitem)
            #mitem = Gtk.SeparatorMenuItem()
            #self.tree_view_menu.append(mitem)

        self.tree_view_menu.show_all()

    
    def open_menu (self, vobject = None, system_id = None):
        """ Function doc """
        self.system_id = system_id
        #print (vobject)
        self.tree_view_menu.popup(None, None, None, None, 0, 0)


