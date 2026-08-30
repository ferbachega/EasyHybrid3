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

from util.debug import dprint
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


            #                       TARGET SURFACE OBJECT SELECTOR
            #'''--------------------------------------------------------------------------------------------'''
            # [EN] Added so that the wireframe/opacity/smooth-shading (and,
            # from now on, MEP colour-scale) controls below apply to ONE
            # specific surface VismolObject, chosen here -- instead of every
            # surface object in the session at once, which is what all
            # three handlers used to do (walking every vm_objects_dic entry
            # with is_surface == True). Repopulated by
            # _refresh_surface_target_combo() every time this window opens
            # and every time a new surface finishes generating (auto-
            # selecting the one just created, so the common "generate then
            # immediately tweak it" flow needs no extra clicking).
            '''
            self.label_surface_target = Gtk.Label(label="Target surface:")
            self.cbx_surface_target = Gtk.ComboBoxText()
            self.cbx_surface_target.set_tooltip_text(
                "Which already-generated surface the controls below (Wireframe, "
                "Opacity, Smooth shading, and the MEP color scale) affect. "
                "Only this one, not every surface in the session.")
            self.box_surface_type.pack_start(self.label_surface_target, False, False, 0)
            self.box_surface_type.pack_start(self.cbx_surface_target, False, False, 0)
            '''
            #'''--------------------------------------------------------------------------------------------'''


            #                       RENDER MODE (filled surface vs wireframe)
            #'''--------------------------------------------------------------------------------------------'''
            self.chk_surface_wireframe = Gtk.CheckButton(label="Wireframe")
            self.chk_surface_wireframe.connect("toggled", self.on_surface_wireframe_toggled)
            #self.box_surface_type.pack_start(self.chk_surface_wireframe, False, False, 0)
            #'''--------------------------------------------------------------------------------------------'''


            #                       OPACITY (surface transparency)
            #'''--------------------------------------------------------------------------------------------'''
            self.label_surface_opacity = Gtk.Label(label="Opacity:")
            self.scale_surface_opacity = Gtk.Scale.new_with_range ( Gtk.Orientation.HORIZONTAL, 0, 100, 1 )
            self.scale_surface_opacity.set_value ( 100 )   # 100% = opaco, igual ao comportamento de antes
            self.scale_surface_opacity.set_size_request ( 120, -1 )
            self.scale_surface_opacity.set_digits ( 0 )
            self.scale_surface_opacity.set_value_pos ( Gtk.PositionType.RIGHT )
            self.scale_surface_opacity.connect ( "value-changed", self.on_surface_opacity_changed )
            self.scale_surface_opacity.set_tooltip_text (
                "Opacity of the surfaces (100% = opaque, 0% = fully "
                "transparent). Applies to every surface already created "
                "in this session (orbitals, potential, density, MEP...).")

            #self.box_surface_type.pack_start(self.label_surface_opacity, False, False, 0)
            #self.box_surface_type.pack_start(self.scale_surface_opacity, True, True, 0)
            #'''--------------------------------------------------------------------------------------------'''


            #                       SHADING (flat vs smooth normals)
            #'''--------------------------------------------------------------------------------------------'''
            self.chk_surface_smooth = Gtk.CheckButton(label="Smooth shading")
            self.chk_surface_smooth.connect("toggled", self.on_surface_smooth_toggled)
            self.chk_surface_smooth.set_tooltip_text(
                "Off (default): constant normal per face (flat "
                "shading) -- shows the marching-cubes triangulation's "
                "facets.\n"
                "On: normal interpolated per vertex, averaged from "
                "adjacent faces (smooth shading) -- surface has a "
                "smoother appearance, with no visible facets.")
            #self.box_surface_type.pack_start(self.chk_surface_smooth, False, False, 0)
            #'''--------------------------------------------------------------------------------------------'''


            #                       MEP COLOR SCALE (vmin/vmax manuais, opcional)
            #'''--------------------------------------------------------------------------------------------'''
            self.label_mep_vmin = Gtk.Label(label="MEP vmin:")
            self.entry_mep_vmin = Gtk.Entry()
            self.entry_mep_vmin.set_width_chars(8)
            self.entry_mep_vmin.set_placeholder_text("auto")
            self.entry_mep_vmin.set_tooltip_text(
                "Minimum electrostatic potential value (in atomic "
                "units, hartree/e) used for the MEP color scale.\n"
                "Empty = automatic (2nd percentile of the computed values).")

            self.label_mep_vmax = Gtk.Label(label="vmax:")
            self.entry_mep_vmax = Gtk.Entry()
            self.entry_mep_vmax.set_width_chars(8)
            self.entry_mep_vmax.set_placeholder_text("auto")
            self.entry_mep_vmax.set_tooltip_text(
                "Maximum electrostatic potential value (in atomic "
                "units, hartree/e) used for the MEP color scale.\n"
                "Empty = automatic (98th percentile of the computed values).")

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
                "Colormap used to map the electrostatic potential to color "
                "(see COLOR_MAPS in util/colormaps.py). Diverging maps "
                "(coolwarm, vik, berlin, bam) follow the red="
                "negative/blue=positive convention; sequential maps (jet, "
                "rainbow, gnuplot, magma, viridis...) have no defined center.")

            self._mep_cmap_names = list ( COLOR_MAPS.keys() )
            for i, name in enumerate ( self._mep_cmap_names ):
                self.cbx_mep_cmap.insert ( i, str(i), name )
            try:
                default_idx = self._mep_cmap_names.index ( 'jet' )
            except ValueError:
                default_idx = 0
            self.cbx_mep_cmap.set_active ( default_idx )

            self.box_surface_type.pack_start(self.label_mep_cmap, False, False, 0)
            self.box_surface_type.pack_start(self.cbx_mep_cmap, False, False, 0)
            self.label_mep_cmap.hide()
            self.cbx_mep_cmap.hide()
            #'''--------------------------------------------------------------------------------------------'''


            #                       MEP POTENTIAL GRID SPACING (optional -- decouples from the density spacing)
            #'''--------------------------------------------------------------------------------------------'''
            # [EN] See changelog item 14: generator.GridPotential() was
            # measured (live, in the user's environment) to take ~180x
            # longer than generator.GridDensity() on the SAME grid --
            # confirmed, by reading pDynamo3's own QCModelBase.py, to be a
            # real, inherent cost difference (GridPotential needs an
            # attraction-integral-like evaluation between every PAIR of
            # basis functions at each grid point, O(n_basis^2) per point,
            # vs GridDensity's O(n_basis) per point). Since MEP now does
            # real trilinear interpolation (_reconstruct_regular_grid +
            # _trilinear_interpolate), the potential does NOT need to be
            # evaluated on the same fine grid as the density surface's own
            # geometry -- a much coarser, cheaper potential grid still
            # interpolates smoothly onto the finer density mesh.
            self.label_mep_pot_spacing = Gtk.Label(label="Potential spacing:")
            self.entry_mep_pot_spacing = Gtk.Entry()
            self.entry_mep_pot_spacing.set_width_chars(8)
            self.entry_mep_pot_spacing.set_placeholder_text("auto (2.5x)")
            self.entry_mep_pot_spacing.set_tooltip_text(
                "Grid spacing used ONLY to calculate the electrostatic "
                "potential (in Bohr -- atomic units, same convention "
                "as the main spacing field). Can be considerably "
                "coarser than the density spacing (which defines the "
                "mesh geometry) with no perceptible visual loss, "
                "because the value is trilinearly interpolated at the "
                "density mesh vertices anyway.\n"
                "Empty = automatic (2.5x the main spacing).\n"
                "Reason: computing the electrostatic potential is MUCH "
                "more expensive than computing the density on the same "
                "grid (measured: ~180x slower in a real case) -- a "
                "coarser grid just for the potential drastically "
                "reduces this cost (2.5x coarser = grid with ~2.5^3 "
                "= ~15x fewer points = ~15x faster).")
            self.box_surface_type.pack_start(self.label_mep_pot_spacing, False, False, 0)
            self.box_surface_type.pack_start(self.entry_mep_pot_spacing, False, False, 0)
            self.label_mep_pot_spacing.hide()
            self.entry_mep_pot_spacing.hide()
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
                dprint(col)
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
            
            # [EN] "Export orbital energies..." button -- user request:
            # print/export the orbital energies (already computed in
            # generate_wavefunction_parallel() below, via
            # system.scratch.orbitalsP.energies) as a plain tab-separated
            # text table, pasteable directly into a spreadsheet program.
            # Added programmatically (no .glade edit) -- same pattern
            # already used throughout this window for the MEP fields,
            # the "Target surface" combobox, etc. selection_treeview (the
            # widget just built above) sits inside an UNNAMED
            # GtkScrolledWindow in the .glade file, itself a sibling of
            # the box holding box_surface_type/btn_external_file/etc --
            # navigated to here via get_parent() twice rather than adding
            # a new named container to the .glade file.
            self.btn_export_orbitals = Gtk.Button ( label = "Export orbital energies..." )
            self.btn_export_orbitals.set_tooltip_text (
                "Saves a table (TAB-separated text, pastes directly into a "
                "spreadsheet) with index, label (HOMO/LUMO), occupancy "
                "and energy (Hartree and eV) for each orbital, for every "
                "frame already computed (every frame with an imported "
                "wavefunction, not just the current frame)." )
            self.btn_export_orbitals.connect ( "clicked", self.on_button_export_orbitals )
            _scrolled_parent = self.treeview.get_parent ( )
            if _scrolled_parent is not None:
                _outer_box = _scrolled_parent.get_parent ( )
                if _outer_box is not None and hasattr ( _outer_box, "pack_start" ):
                    _outer_box.pack_start ( self.btn_export_orbitals, False, False, 4 )
                    self.btn_export_orbitals.show ( )
            # [EN] User request: "ao clicar em um orbital da treeview, ele
            # printe os dados da LCAO no terminal, por hora." --
            # "cursor-changed" fires on a plain single click/selection
            # change (unlike "row-activated", which normally needs a
            # double-click or Enter) -- matches "ao clicar" literally.
            self.treeview.connect ( "cursor-changed", self.on_orbital_row_selected )
            #-----------------------------------------------------------
            self.btn_render = self.builder.get_object('btn_render')
            self.btn_render.connect('clicked', self.on_render_button)
            #-----------------------------------------------------------
            
            
            self.btn_color_plus  = self.builder.get_object('btn_color_plus')
            self.btn_color_minus = self.builder.get_object('btn_color_minus')
            #self.btn_color_density = self.builder.get_object('btn_color_density')
            # Setting a specific color (for example, red)
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
                "Density or orbital .cube file (e.g. generated by ORCA's "
                "orca_plot). The surface mesh comes from here.")

            self.label_external_potential = Gtk.Label(label="Potential .cube (optional, for MEP):")
            #self.btn_external_potential_file = Gtk.FileChooserButton(
            #    title  = "Select the potential .cube file (optional)",
            #    action = Gtk.FileChooserAction.OPEN,
            #)
            _cube_filter = Gtk.FileFilter()
            _cube_filter.set_name("Cube files (*.cube)")
            _cube_filter.add_pattern("*.cube")
            #self.btn_external_potential_file.add_filter(_cube_filter)
            #self.btn_external_potential_file.connect("file-set", self.on_external_potential_file_set)
            #self.btn_external_potential_file.set_tooltip_text(
            #    "Optional. If given, colors the density/orbital mesh "
            #    "above by electrostatic potential interpolated from this "
            #    "cube (a continuous color map, same as MEP -- see mep_colormap).")
            self.box_surface_type.pack_start(self.label_external_potential, False, False, 0)
            #self.box_surface_type.pack_start(self.btn_external_potential_file, False, False, 0)
            self.label_external_potential.hide()
            #self.btn_external_potential_file.hide()



            self.builder.get_object('orbital_scrolled_window').set_size_request(-1, 300)

            #'''--------------------------------------------------------------------------------------------'''

            # so agora, com TODOS os widgets acima ja criados, e que ativamos o
            # item 0 do combobox -- set_active() dispara "changed" (chama
            # surface_combobox_change) IMMEDIATELY and synchronously, so it
            # must come after everything this handler may try to
            # show/hide (otherwise AttributeError: widget does not exist yet).
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
    #
    # [EN] UPDATED: these three handlers used to act on EVERY surface
    # object in the session at once -- the user reported this directly
    # ("quando optamos entre lines ou triangulos... altera a representacao
    # de TODAS as superficies"). They now act on exactly ONE object: the
    # one currently chosen in self.cbx_surface_target (see that combobox's
    # own comment, and _refresh_surface_target_combo() /
    # _get_target_surface_object() right below).
    def _refresh_surface_target_combo (self, select_index = None):
        """ Repopulates self.cbx_surface_target from every vismol object in
        the session currently flagged is_surface == True. Call this: (a)
        once when the window opens, so surfaces from an earlier session
        already show up; (b) right after a new surface finishes
        generating, passing select_index=<the new object's own .index>,
        so it becomes the active target immediately (the common
        "generate then tweak it" flow needs no extra clicking to select
        what was just made). """
        self.cbx_surface_target.remove_all ( )
        surfaces = [ v for v in self.vm_session.vm_objects_dic.values ( )
                     if getattr ( v, "is_surface", False ) ]
        active_row = 0
        for row, vobject in enumerate ( surfaces ):
            self.cbx_surface_target.append ( str ( vobject.index ), vobject.name )
            if select_index is not None and vobject.index == select_index:
                active_row = row
        if surfaces:
            self.cbx_surface_target.set_active ( active_row )

    def _get_target_surface_object (self):
        """ Returns the VismolObject currently chosen in
        self.cbx_surface_target, or None if there isn't one (combobox
        empty -- no surface generated yet in this session). """
        active_id = self.cbx_surface_target.get_active_id ( )
        if active_id is None:
            return None
        return self.vm_session.vm_objects_dic.get ( int ( active_id ) )

    def on_surface_wireframe_toggled (self, widget):
        """ Alterna entre superficie preenchida (GL_FILL) e wireframe (GL_LINE)
        para a superficie escolhida em self.cbx_surface_target (so essa --
        ver o comentario acima e o da propria combobox). Nao recalcula
        malha nem recompila shader -- so muda o render_mode de cada
        representacao encontrada nesse UM objeto (ver representations.py). """
        mode = "lines" if widget.get_active() else "surface"
        vobject = self._get_target_surface_object ( )
        if vobject is None:
            return
        for rep in vobject.representations.values():
            if hasattr(rep, "set_render_mode"):
                rep.set_render_mode(mode)
        self.vm_session.vm_glcore.queue_draw()

    def on_surface_opacity_changed (self, widget):
        """ Ajusta a opacidade (alpha) da superficie escolhida em
        self.cbx_surface_target (so essa). widget.get_value() vai de 0
        a 100 (%); SurfaceRepresentation.set_alpha() espera 0.0-1.0. """
        alpha = widget.get_value() / 100.0
        vobject = self._get_target_surface_object ( )
        if vobject is None:
            return
        for rep in vobject.representations.values():
            if hasattr(rep, "set_alpha"):
                rep.set_alpha(alpha)
        self.vm_session.vm_glcore.queue_draw()

    def on_surface_smooth_toggled (self, widget):
        """ Alterna entre flat shading (normal por face) e smooth shading
        (normal por vertice, media das faces adjacentes) na superficie
        escolhida em self.cbx_surface_target (so essa). Nao recalcula a
        malha -- as normais suaves ja foram calculadas na geracao
        (surface.MakeVertexNormalsFromPolygonalNormals(), nativo do
        pDynamo3) e enviadas pro VAO; aqui so muda qual delas o shader
        usa (ver geometry_shader_surface). """
        mode = "smooth" if widget.get_active() else "flat"
        vobject = self._get_target_surface_object ( )
        if vobject is None:
            return
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
        dprint(index)

        if index == 4:
            self.label_mep_vmin.show()
            self.entry_mep_vmin.show()
            self.label_mep_vmax.show()
            self.entry_mep_vmax.show()
            self.label_mep_cmap.show()
            self.cbx_mep_cmap.show()
            self.label_mep_pot_spacing.show()
            self.entry_mep_pot_spacing.show()
        else:
            self.label_mep_vmin.hide()
            self.entry_mep_vmin.hide()
            self.label_mep_vmax.hide()
            self.entry_mep_vmax.hide()
            self.label_mep_cmap.hide()
            self.cbx_mep_cmap.hide()
            self.label_mep_pot_spacing.hide()
            self.entry_mep_pot_spacing.hide()

        if index == 1:
            #self.builder.get_object('btn_import_wavefunction').set_sensitive(False)
            #self.builder.get_object('selection_treeview')     .set_sensitive(False)
            self.builder.get_object('btn_import_wavefunction').hide()
            self.builder.get_object('selection_treeview')     .hide()
            self.builder.get_object('orbital_scrolled_window')     .hide()
            self.builder.get_object('label_external_file').hide()
            self.builder.get_object('btn_external_file')  .hide()
            self.label_external_potential.hide()
            
            self.builder.get_object('btn_color_minus').show()
            self.builder.get_object('btn_color_plus' ).show()
            self.builder.get_object('label_color_minus' ).show()
            self.builder.get_object('label_color_plus' ).show()
            self.btn_export_orbitals.hide()
            #self.btn_external_potential_file.hide()
            self.window.queue_resize()
            self.window.resize(1,1)
            
        elif index == 2:
            #self.builder.get_object('btn_import_wavefunction').set_sensitive(False)
            #self.builder.get_object('selection_treeview')     .set_sensitive(False)
            self.builder.get_object('btn_import_wavefunction').hide()
            self.builder.get_object('selection_treeview')     .hide()
            self.builder.get_object('orbital_scrolled_window')     .hide()
            self.builder.get_object('label_external_file').hide()
            self.builder.get_object('btn_external_file')  .hide()
            self.label_external_potential.hide()
            
            self.builder.get_object('btn_color_minus').hide()
            self.builder.get_object('btn_color_plus' ).show()
            self.builder.get_object('label_color_minus' ).hide()
            self.builder.get_object('label_color_plus' ).show()
            self.btn_export_orbitals.hide()
            self.window.queue_resize()
            self.window.resize(1,1)            
        elif index == 3:
            self.builder.get_object('label_external_file').show()
            self.builder.get_object('btn_external_file').show()
            self.label_external_potential.show()
            #self.btn_external_potential_file.show()
            
            #self.builder.get_object('btn_import_wavefunction').set_sensitive(False)
            #self.builder.get_object('orbital_scrolled_window'     ).set_sensitive(False)
            #self.builder.get_object('selection_treeview'     ).set_sensitive(False)
            
            self.builder.get_object('btn_import_wavefunction').hide()
            self.builder.get_object('selection_treeview'     ).hide()
            
            self.builder.get_object('btn_color_minus').show()
            self.builder.get_object('btn_color_plus' ).show()
            self.builder.get_object('label_color_minus' ).show()
            self.builder.get_object('label_color_plus' ).show()
            self.builder.get_object('btn_import_wavefunction').hide()
            self.builder.get_object('selection_treeview')     .hide()
            
            self.btn_export_orbitals.hide()
            self.window.queue_resize()
            self.window.resize(1,1)            
        elif index == 4:
            #self.builder.get_object('btn_import_wavefunction').set_sensitive(False)
            #self.builder.get_object('selection_treeview')     .set_sensitive(False)
            self.builder.get_object('btn_import_wavefunction').hide()
            self.builder.get_object('selection_treeview')     .hide()
            self.builder.get_object('orbital_scrolled_window')     .hide()
            self.builder.get_object('label_external_file').hide()
            self.builder.get_object('btn_external_file')  .hide()
            self.label_external_potential.hide()
            
            self.builder.get_object('btn_color_minus').hide()
            self.builder.get_object('btn_color_plus' ).hide()
            self.builder.get_object('label_color_minus' ).hide()
            self.builder.get_object('label_color_plus' ).hide()
            
            self.btn_export_orbitals.hide()
            self.window.queue_resize()
            self.window.resize(1,1)
        else:
            self.builder.get_object('btn_import_wavefunction').set_sensitive(True)
            self.builder.get_object('selection_treeview')     .set_sensitive(True)
            self.builder.get_object('btn_import_wavefunction').show()
            self.builder.get_object('orbital_scrolled_window')     .show()
            self.builder.get_object('selection_treeview')     .show()
            self.builder.get_object('label_external_file').hide()
            self.builder.get_object('btn_external_file')  .hide()
            self.label_external_potential.hide()

            self.builder.get_object('btn_color_minus').show()
            self.builder.get_object('btn_color_plus' ).show()
            self.builder.get_object('label_color_minus' ).show()
            self.builder.get_object('label_color_plus' ).show()
            
            self.btn_export_orbitals.show()
            self.window.queue_resize()
            self.window.resize(1,1) 

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
                dprint('vobject has no Normal Modes data')
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
        dprint(key, vismol_object.frames.shape[0], model, model[iter][1])
        name = str(key) +' '+model[iter][1]#+' '+ str(model[iter][3])
        #_GridSpacing = 0.6
        _OrbitalTag    = "Grid Orbitals"
        _IsosurfaceTag = "Isosurface"
        
        
        trajectory = [None]*vismol_object.frames.shape[0]
        joblist = []
        
        for frame in range(vismol_object.frames.shape[0]):
            #self.p_session.set_psystem_coordinates_from_vobject(vobject)
            #'''
            dprint(vismol_object, frame)
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
        vobject_tmp.surface_type = "orbital"   # usado por _surf_setup() no menu de contexto da arvore
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
        
        dprint (results)
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
        dprint(vobject_id)
        
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
        #rgba_density = self.btn_color_density.get_rgba()
        rgba_density = self.btn_color_plus.get_rgba()
        
        color_density =  [rgba_density.red , rgba_density.green,  rgba_density.blue ] 
        color_plus  = [rgba_plus.red , rgba_plus.green,  rgba_plus.blue ]
        color_minus = [rgba_minus.red, rgba_minus.green, rgba_minus.blue]
        
        index = self.cbx_surface_type.get_active()
        #print(index, color_minus, color_plus)

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
                'color_plus'     : color_density ,
                'color_minus'    : color_density ,
                
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
                                                                                iso_color     = color_density             ,
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
            vobject_tmp.surface_type = "density"   # usado por _surf_setup() no menu de contexto da arvore
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            dprint('\n\nvismol_object.e_treeview_iter', vismol_object.e_treeview_iter,'\n\n')
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
            _mep_pot_spacing = _parse_optional_float ( self.entry_mep_pot_spacing )
            # None = automatico (2.5x o espacamento principal, ver
            # generate_grid_parallel() branch 'mep' e o changelog item 14).

            joblist = []
            for frame in range(vismol_object.frames.shape[0]):
                self.p_session.set_psystem_coordinates_from_vobject( vobject = vismol_object, 
                                                                               system_id = None, 
                                                                               frame = frame)
                parameters = {
                'type'              : 'mep',
                '_GridSpacing'      : _GridSpacing,
                '_OrbitalTag'       : 'density_mep',
                '_isovalue'         : _isovalue,
                '_IsosurfaceTag'    : 'Isosurface',
                'orbital_key'       : 0,
                'color_plus'        : color_plus  ,
                'color_minus'       : color_minus ,
                'mep_vmin'          : _mep_vmin   ,
                'mep_vmax'          : _mep_vmax   ,
                'mep_cmap_name'     : _mep_cmap_name ,
                'mep_pot_spacing'   : _mep_pot_spacing ,
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
            # mep_colormap in generate_grid_parallel), not used by
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
            vobject_tmp.surface_type = "mep"   # usado por _surf_setup() no menu de contexto da arvore
            # [EN] Caches, per frame, exactly what's needed to redo the
            # colour scale (vmin/vmax/colormap) cheaply later -- see the
            # matching comment where these two keys are added to the
            # returned dict inside generate_grid_parallel()'s 'mep'
            # branch. Kept as one list per frame (mirroring
            # surface_trajectory's own per-frame structure) rather than
            # just the first frame, so a future per-frame MEP recolor
            # (not implemented yet -- current _surf_setup only recolors
            # whichever frame is currently displayed) has what it needs
            # without another full recompute.
            vobject_tmp.mep_potential_values = [ r.get ( 'mep_raw_potential_values' ) for r in results ]
            vobject_tmp.mep_polygons         = [ r.get ( 'mep_raw_polygons' )         for r in results ]
            vobject_tmp.mep_vmin      = _mep_vmin
            vobject_tmp.mep_vmax      = _mep_vmax
            vobject_tmp.mep_cmap_name = _mep_cmap_name
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
            self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
            self.main.refresh_widgets()
            self.vm_session.vm_glcore.queue_draw()
            self.counter +=1

        elif index == 3:
            if not self.external_density_path:
                dprint("No density/orbital .cube file selected -- "
                      "click 'Choose file...' before rendering.")
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

            # an external cube is a static file -- no system/coords/QC
            # involved at all, just file reading + marching cubes.
            # Roda direto (sem multiprocessing.Pool: e so I/O + um
            # compiled algorithm, there is no heavy QC calculation to parallelize
            # aqui, e cada mudanca de aba abriria um Pool novo a toa).
            try:
                single_result = generate_grid_parallel ( [ 0, None, None, parameters ] )
            except CubeFileError as error:
                dprint ( "Error reading .cube file: {}".format ( error ) )
                return False

            # replica o mesmo resultado pra todos os "frames" do objeto
            # parent, just for safety (surface_trajectory[frame] cannot
            # dar index error se o usuario tiver uma trajetoria carregada
            # e trocar de frame -- o cubo externo e sempre a mesma
            # static mesh, does not change per frame).
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
            # [EN] External cube surfaces are only MEP-coloured when a
            # potential .cube was ALSO given (self.external_potential_path);
            # otherwise it is a flat-colour surface like density/orbital.
            # _generate_external_cube_surface() only includes
            # mep_raw_potential_values/mep_raw_polygons in its returned dict
            # for the MEP case (see that function), so results[0].get(...) is
            # None in the non-MEP case and the cache attributes below are
            # simply left unset (getattr(..., None) elsewhere handles that).
            if self.external_potential_path:
                vobject_tmp.surface_type = "mep"   # usado por _surf_setup() no menu de contexto da arvore
                vobject_tmp.mep_potential_values = [ r.get ( 'mep_raw_potential_values' ) for r in results ]
                vobject_tmp.mep_polygons         = [ r.get ( 'mep_raw_polygons' )         for r in results ]
                vobject_tmp.mep_vmin      = parameters.get ( 'mep_vmin' )   # sempre None por enquanto (ver comentario na montagem de parameters, acima)
                vobject_tmp.mep_vmax      = parameters.get ( 'mep_vmax' )
                vobject_tmp.mep_cmap_name = parameters.get ( 'mep_cmap_name', 'coolwarm' )
            else:
                vobject_tmp.surface_type = "external"   # usado por _surf_setup() no menu de contexto da arvore
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
            vobject_tmp.surface_type = "potential"   # usado por _surf_setup() no menu de contexto da arvore
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
            dprint(key, vismol_object.frames.shape[0], model, model[iter][1])
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
            vobject_tmp.surface_type = "orbital"   # usado por _surf_setup() no menu de contexto da arvore
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
        dprint('on_button_import_wavefunction')
        
        system_id = self.system_names_combo.get_system_id()
        system    = self.main.p_session.psystem[system_id]
        
        vobject_id    = self.coordinates_combobox.get_vobject_id()
        vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]
        dprint(system_id, vobject_id, system, vismol_object)
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
        
        dprint()
        self.treeview.set_model(self.orbital_liststore_dict[vobject_id][self.frame])

    def on_orbital_row_selected (self, treeview):
        """ [EN] User request: capture and print the LCAO (linear
        combination of atomic orbitals) coefficients for whichever
        orbital was just clicked in the orbital list -- "por hora"
        (for now) just prints to the console/terminal; a proper GUI
        display is a possible later step.

        Reads orbitals_matrix (the (nBasis, nOrbitals) coefficient
        matrix -- see generate_wavefunction_parallel()'s own comment for
        exactly what this is and how it was captured), center_function_
        pointers (basis-function index boundaries per QC atom), and
        atom_symbols, all cached per-frame in self.wave_function_dict
        alongside the orbital energies themselves.

        Groups the raw coefficients by QC ATOM (not down to individual
        shells/angular-momentum components -- that finer level is
        available in principle, see the comment in
        generate_wavefunction_parallel(), but not implemented here to
        keep this first pass simple) and also prints a rough per-atom
        "weight" (sum of squared coefficients for that atom's basis
        functions) -- explicitly NOT a real Mulliken population (that
        needs the overlap matrix too, to account for the atomic-orbital
        basis not being orthonormal), just a quick, everyday indicator
        of which atoms contribute the most to this particular orbital. """
        selection = treeview.get_selection ( )
        model, treeiter = selection.get_selected ( )
        if treeiter is None:
            return

        orbital_index = model[treeiter][0]   # coluna 0 do liststore = indice do orbital (ver generate_wavefunction_parallel)
        orbital_label  = model[treeiter][1]
        orbital_energy = model[treeiter][3]

        vobject_id = self.coordinates_combobox.get_vobject_id ( )
        wave_data  = self.wave_function_dict.get ( vobject_id )
        if not wave_data:
            dprint ( "No orbital computed yet -- click 'Import Wavefunction' first." )
            return

        frame = self.frame if self.frame < len ( wave_data ) else 0
        frame_data = wave_data[frame]
        if len ( frame_data ) < 6:
            dprint ( "LCAO data not available (computed with an older version, before this capture existed)." )
            return

        orbitals_matrix, center_function_pointers, atom_symbols = frame_data[3], frame_data[4], frame_data[5]
        if orbitals_matrix is None or center_function_pointers is None or atom_symbols is None:
            print ( "LCAO data not available for this frame (see the error message from the wavefunction import)." )
            return

        coeffs = orbitals_matrix[:, orbital_index]

        print ( )
        print ( "=== LCAO do orbital {} ({}, frame {}, energia = {:.6f} Hartree) ===".format (
                 orbital_index, orbital_label, frame, orbital_energy ) )
        for atom_idx, symbol in enumerate ( atom_symbols ):
            start = int ( center_function_pointers[atom_idx] )
            end   = int ( center_function_pointers[atom_idx + 1] )
            atom_coeffs = coeffs[start:end]
            rough_weight = float ( np.sum ( atom_coeffs ** 2 ) )
            coeffs_str = ", ".join ( "{:+.4f}".format ( c ) for c in atom_coeffs )
            print ( "  Atom {:3d} ({:>2s})  weight~{:6.2%}  coefficients = [{}]".format (
                     atom_idx, symbol, rough_weight, coeffs_str ) )
        print ( "  (weight = sum of squared coefficients per atom -- this is NOT a real "
                 "Mulliken population, which would also need the overlap matrix; "
                 "it is just a quick indicator of which atoms contribute the most.)" )

    def on_button_export_orbitals (self, widget):
        """ [EN] Writes self.wave_function_dict[vobject_id] (the orbital
        energies already computed by on_button_import_wavefunction() /
        generate_wavefunction_parallel() above -- system.scratch.
        orbitalsP.energies, one list of [index, label, occupancy, energy,
        False] per frame) out as a plain TAB-separated text table, for
        EVERY frame that has been computed (not just the one currently
        shown in self.treeview) -- pastes directly into a spreadsheet
        program (Excel, LibreOffice Calc, Google Sheets: TAB-separated
        text pastes as separate columns in all of them, no import wizard
        needed).

        Energy is given in both Hartree (pDynamo's native unit -- see
        orbitalsP.energies itself, never converted anywhere else in this
        file) and eV (the unit HOMO/LUMO gaps are conventionally
        discussed in), so the exported table is directly usable without
        the user needing to convert by hand.

        Layout: each frame gets its OWN block of 5 columns (Index,
        Label, Occupancy, Energy Hartree, Energy eV), placed SIDE BY
        SIDE with the other frames' blocks (not stacked one below the
        other) -- makes it easy to compare the same orbital across
        frames just by reading across a row in the spreadsheet, rather
        than scrolling down through one long list per frame. Rows
        within each block are in natural (ascending) orbital-index
        order -- NOT the reversed order self.treeview displays on
        screen (highest index first). """
        vobject_id = self.coordinates_combobox.get_vobject_id ( )
        wave_data  = self.wave_function_dict.get ( vobject_id )
        if not wave_data:
            dprint ( "No orbital computed yet for this object -- "
                     "click 'Import Wavefunction' first." )
            return

        HARTREE_TO_EV = 27.211386245988   # constante fisica padrao (CODATA)

        # [EN] Layout changed at the user's request: was one row per
        # (frame, orbital) pair, stacking every frame's block UNDERNEATH
        # the previous one (long, hard to compare frames side by side).
        # Now each frame gets its own group of 5 columns (Index, Label,
        # Occupancy, Energy Hartree, Energy eV), placed NEXT TO the
        # previous frame's block instead -- one "Frame N" header row
        # above each group, then a shared sub-header row, then one data
        # row per orbital index. Assumes every frame has the same number
        # of orbitals in the same order (true for the same QC system/
        # basis set across frames of one trajectory) -- uses the
        # SMALLEST orbital count found, defensively, in case that's ever
        # not the case, rather than raising an index error.
        n_frames   = len ( wave_data )
        n_orbitals = min ( len ( data[0] ) for data in wave_data ) if wave_data else 0

        header_frame_row = []
        header_column_row = []
        for frame in range ( n_frames ):
            header_frame_row.append ( "Frame {}".format ( frame ) )
            header_frame_row.extend ( [ "" ] * 4 )   # preenche o resto do bloco de 5 colunas, so pro rotulo "Frame N" ficar alinhado com a primeira coluna do bloco
            header_column_row.extend ( [ "Index", "Label", "Occupancy", "Energy (Hartree)", "Energy (eV)" ] )

        lines = [ "\t".join ( header_frame_row ), "\t".join ( header_column_row ) ]

        for i in range ( n_orbitals ):
            row_parts = []
            for frame in range ( n_frames ):
                orbitals = wave_data[frame][0]
                idx, label, occupancy, energy, _flag = orbitals[i]
                energy_ev = energy * HARTREE_TO_EV
                row_parts.extend ( [ str ( idx ), str ( label ), str ( occupancy ),
                                      "{:.8f}".format ( energy ), "{:.6f}".format ( energy_ev ) ] )
            lines.append ( "\t".join ( row_parts ) )

        # [EN] User request: "armazene essa informacao nos dados dos
        # orbitais (pode ficar no txt gerado tambem)" -- adds a SECOND
        # table below the energy one, one row per orbital, giving each QC
        # atom's rough LCAO weight (sum of squared coefficients for that
        # atom's basis functions -- see on_orbital_row_selected()'s own
        # comment for why this is a quick indicator and NOT a true
        # Mulliken population). One column per atom (not per individual
        # basis function/coefficient -- that would make the table
        # impractically wide for anything but the smallest molecules) --
        # kept to the FIRST frame that actually has LCAO data captured
        # (atom count/order is the same across every frame of one
        # trajectory anyway, so repeating it per frame would just be
        # redundant column groups).
        lcao_frame = None
        for frame_data in wave_data:
            if len ( frame_data ) >= 6 and frame_data[3] is not None:
                lcao_frame = frame_data
                break

        if lcao_frame is not None:
            orbitals_matrix, center_function_pointers, atom_symbols = lcao_frame[3], lcao_frame[4], lcao_frame[5]
            lines.append ( "" )   # linha em branco separando as duas tabelas
            lines.append ( "LCAO atom weights (sum of squared coefficients per atom -- "
                            "NOT a real Mulliken population, see on_orbital_row_selected)" )
            atom_header = [ "Index", "Label" ] + [ "Atom{} ({})".format ( a, s ) for a, s in enumerate ( atom_symbols ) ]
            lines.append ( "\t".join ( atom_header ) )
            for i in range ( n_orbitals ):
                idx, label = wave_data[0][0][i][0], wave_data[0][0][i][1]
                coeffs = orbitals_matrix[:, i]
                row_parts = [ str ( idx ), str ( label ) ]
                for atom_idx in range ( len ( atom_symbols ) ):
                    start = int ( center_function_pointers[atom_idx] )
                    end   = int ( center_function_pointers[atom_idx + 1] )
                    weight = float ( np.sum ( coeffs[start:end] ** 2 ) )
                    row_parts.append ( "{:.6f}".format ( weight ) )
                lines.append ( "\t".join ( row_parts ) )

        text_content = "\n".join ( lines )

        dialog = Gtk.FileChooserDialog (
            title  = "Export orbital energies",
            parent = self.window,
            action = Gtk.FileChooserAction.SAVE,
        )
        dialog.add_buttons (
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_SAVE,   Gtk.ResponseType.OK,
        )
        dialog.set_do_overwrite_confirmation ( True )
        dialog.set_current_name ( "orbital_energies.txt" )
        response = dialog.run ( )
        if response == Gtk.ResponseType.OK:
            filepath = dialog.get_filename ( )
            if not filepath.lower ( ).endswith ( ".txt" ):
                filepath += ".txt"
            with open ( filepath, "w" ) as f:
                f.write ( text_content )
            dprint ( "Orbital energies exported to:", filepath )
            dialog.destroy ( )
            # [EN] BUG FIX (user asked "onde ele salvou o
            # orbital_energies.txt?" -- the only confirmation was the
            # print() above, invisible if EasyHybrid was launched from a
            # desktop/menu shortcut instead of a terminal). Shows the
            # exact saved path in an actual dialog too, so it's visible
            # regardless of how the app was launched.
            confirm = Gtk.MessageDialog (
                transient_for = self.window,
                flags         = 0,
                message_type  = Gtk.MessageType.INFO,
                buttons       = Gtk.ButtonsType.OK,
                text          = "Orbital energies exported.",
            )
            confirm.format_secondary_text ( filepath )
            confirm.run ( )
            confirm.destroy ( )
        else:
            dialog.destroy ( )

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
            dprint ('event.button == 1')



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
            "pDynamo's RegularGrid only accepts axis-aligned grids "
            "(diagonal voxel vectors, no rotation/shear). This "
            ".cube has non-orthogonal axes -- uncommon, but some programs "
            "can generate this; not supported by this importer."
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
    # pDynamo, so a simple C-order reshape is enough.
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


def _pdynamo_array_to_numpy ( arr, dtype ):
    """ [EN] Converts a pDynamo3 Array2D-like object (polygons, vertices,
    vertexNormals, ...) to a plain numpy array as fast as possible, so the
    rest of the pipeline can use vectorised numpy operations instead of
    per-element Python loops (see changelog item 13 -- this whole family
    of helper functions exists because per-triangle/per-vertex Python
    loops in surface_parser() / surface_parser_mep() /
    _compute_valid_polygon_mask() turned out to be the single largest
    bottleneck in surface generation, benchmarked at ~15-45x slower than
    the vectorised equivalent on a realistic ~30k-vertex/60k-triangle
    mesh).

    Tries np.asarray() first (works for free if pDynamo3's Array2D
    implements the buffer/array protocol -- COULD NOT BE CONFIRMED in the
    assistant's environment, since no pDynamo3 installation was available
    to test against; this is an optimistic fast path, not a verified
    one). Falls back to a per-ROW (not per-element) Python loop otherwise,
    which is still correct and still meaningfully faster than the
    previous per-element indexing pattern, just not as fast as a native
    numpy conversion would be. Either path produces an identical,
    correct result -- only the speed differs. """
    try:
        result = np.asarray ( arr, dtype = dtype )
        if result.ndim == 2 and result.shape[0] == arr.rows:
            return result
    except Exception:
        pass
    n = arr.rows
    first_row = arr[0, :]
    ncols = len ( first_row )
    out = np.empty ( (n, ncols), dtype = dtype )
    for i in range ( n ):
        out[i, :] = arr[i, :]
    return out


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
    entrada por linha de `polygons` (True = triangulo valido, mantem).

    [EN] VECTORISED (see changelog item 13) -- was a per-triangle Python
    loop before, benchmarked ~44x slower than this on a 60k-triangle mesh. """
    n_tri = polygons.rows
    if n_tri == 0:
        return np.zeros ( 0, dtype = bool ), 0

    polygons_np = _pdynamo_array_to_numpy ( polygons, np.int64 )
    vertices_np = _pdynamo_array_to_numpy ( vertices, np.float64 )
    tri_verts = vertices_np[polygons_np]          # (n_tri, 3, 3) -- fancy indexing, uma leitura so
    e01 = np.linalg.norm ( tri_verts[:,1,:] - tri_verts[:,0,:], axis = 1 )
    e12 = np.linalg.norm ( tri_verts[:,2,:] - tri_verts[:,1,:], axis = 1 )
    e20 = np.linalg.norm ( tri_verts[:,0,:] - tri_verts[:,2,:], axis = 1 )
    max_edge = np.maximum ( np.maximum ( e01, e12 ), e20 )

    typical = np.median ( max_edge )
    if typical == 0.0:
        typical = 1e-9   # fully degenerate mesh -- avoids division by zero later, without breaking
    threshold = size_factor * typical
    mask = max_edge <= threshold
    n_discarded = int ( n_tri - mask.sum ( ) )
    return mask, n_discarded


def surface_parser ( surface, iso_color):
    # [EN] VECTORISED (changelog item 13) -- was a per-triangle, per-vertex,
    # per-component Python loop before (~1000s of tiny list.append() calls),
    # benchmarked ~18x slower than this on a realistic ~60k-triangle mesh.
    # Correctness verified against the old loop-based version on synthetic
    # data before replacing it (identical output, byte for byte).
    """ Function doc """
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
    if skipped:
        dprint ( "surface_parser: {} triangulo(s) fantasma (aresta anormal) descartado(s)".format ( skipped ) )

    polygons_np = _pdynamo_array_to_numpy ( polygons, np.int64 )[valid_mask]
    vertices_np = _pdynamo_array_to_numpy ( vertices, np.float64 )
    normals_np  = _pdynamo_array_to_numpy ( smooth_normals, np.float64 )

    tri_verts = vertices_np[polygons_np] / 1.889725989   # (n_valid_tri, 3, 3) -- Bohr -> Angstrom, fancy indexing (uma leitura vetorizada)
    tri_norms = normals_np[polygons_np]                  # (n_valid_tri, 3, 3)

    vertices_out = tri_verts.reshape ( -1 ).astype ( np.float32 )
    normals_out  = tri_norms.reshape ( -1 ).astype ( np.float32 )

    n_valid_tri = polygons_np.shape[0]
    colors_out  = np.tile ( np.asarray ( iso_color, dtype = np.float32 ), n_valid_tri * 3 )

    # one index per vertex (not per component/float) -- the old version
    # gerava indexes com 3x mais entradas do que vertices de verdade
    # existem no buffer (passava despercebido, mas era um out-of-bounds
    # read em potencial na GPU -- ver nota no README).
    indexes_out = np.arange ( vertices_out.shape[0] // 3, dtype = np.uint32 )
    return vertices_out, colors_out, indexes_out, normals_out


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
    # beyond +-limit (the outliers cut by the percentile) are
    # clamped to t=0 or t=1 by np.clip -- they stay saturated, not
    # quebram a escala.
    t = 0.5 + 0.5 * np.clip ( values / limit, -1.0, 1.0 )

    if reverse:
        t = 1.0 - t

    return _colormap_lookup ( t, COLOR_MAPS[cmap_name] )


# ============================================================================
#  Cheap recolouring, WITHOUT touching pDynamo/marching cubes/multiprocessing
# ============================================================================
# [EN] Both functions below exist to answer the user's original question
# directly ("what information can we store in the vobject so we do not need
# fazer todos o calculo de superficies novamente?"): the mesh geometry
# (vertices/indexes/normals) for a given surface never needs to change
# just because its DISPLAY colour does. VismolObject.surface_trajectory
# (already existing, populated when the surface was first generated --
# see on_render_button() above) stores, per frame, per named lobe
# ("obital_plus"/"obital_minus"), the exact 4-tuple (vertices, colors,
# indexes, normals) that SurfaceRepresentation.draw_representation()
# re-reads on EVERY single frame (unconditionally -- see representations.py).
# That means replacing just the "colors" entry of that cached tuple and
# calling queue_draw() is enough to change what's on screen -- no need to
# touch surfacetrajectory's vertices/indexes/normals, no need to recreate
# the representation, and (for MEP specifically) no need to re-run
# GridDensity/Isosurface/GridPotential/build_potential_interpolator or any
# of the other genuinely expensive steps in generate_grid_parallel()'s
# 'mep' branch (measured there: GridPotential alone can be ~180x the cost
# of GridDensity on the same grid -- see changelog item 14).
#
# Used by _surf_setup() in gui/main/treeview_menu.py.

def recolor_surface_lobe ( vobject, surf_name, new_rgb ):
    """ Cheap recolour for a single-colour-per-lobe surface (orbital,
    density, or potential -- the ones with a "color_plus"/"color_minus"
    pair, one flat colour repeated across every vertex of that lobe).
    Just re-tiles new_rgb across however many vertices that lobe's
    CACHED mesh already has -- doesn't need the vertex count from
    anywhere else, doesn't touch pDynamo, doesn't re-run marching cubes.

    vobject   : the surface VismolObject (vobject.surface_type in
                ("orbital", "density", "potential"), though this function
                itself doesn't check that -- it only needs
                vobject.surface_trajectory[frame][surf_name] to exist).
    surf_name : "obital_plus" or "obital_minus" (matches the key used
                when the surface was first generated -- see
                SurfaceRepresentation(..., surface_name=...) calls in
                on_render_button()).
    new_rgb   : an (r, g, b) sequence, each 0.0-1.0.

    Updates every frame, not just the currently-displayed one (a static,
    single-frame surface is by far the common case, but this stays
    correct even for a multi-frame one). """
    new_rgb = np.asarray ( new_rgb, dtype = np.float32 )
    for frame_data in vobject.surface_trajectory:
        if frame_data is None or surf_name not in frame_data:
            continue
        vertices, colors, indexes, normals = frame_data[surf_name]
        n_vertices = vertices.shape[0] // 3   # vertices e um array 1-D achatado (x,y,z,x,y,z,...)
        new_colors = np.tile ( new_rgb, n_vertices ).astype ( np.float32 )
        frame_data[surf_name] = [ vertices, new_colors, indexes, normals ]


def recolor_mep_surface ( vobject, vmin = None, vmax = None, cmap_name = None ):
    """ Cheap recolour for an MEP surface: redoes ONLY mep_colormap() (the
    scalar-potential -> RGB mapping) plus the per-vertex -> per-triangle-
    corner colour expansion (vertex_colors[polygons]) -- using the RAW,
    pre-colormap potential values and the mask-filtered polygon index
    array, both cached on the vobject the first time it was generated
    (see the 'mep' branch of generate_grid_parallel(), and the External-
    cube-with-MEP branch of _generate_external_cube_surface(), both in
    this same file). Does NOT touch the mesh geometry, does NOT
    re-evaluate the potential at any vertex, does NOT re-run marching
    cubes or any pDynamo QC step.

    vmin, vmax, cmap_name: same meaning as mep_colormap()'s own
    parameters -- pass None to keep vobject's CURRENT value for that one
    (so a dialog can change just vmin while leaving vmax/cmap_name alone,
    for instance). Whatever values end up being used are written back to
    vobject.mep_vmin/mep_vmax/mep_cmap_name, so the setup dialog can show
    the actual current state next time it's opened.

    Raises ValueError if vobject has no cached mep_potential_values/
    mep_polygons (e.g. it's a plain, non-MEP surface, or was generated
    before this caching existed). """
    potential_values_per_frame = getattr ( vobject, "mep_potential_values", None )
    polygons_per_frame         = getattr ( vobject, "mep_polygons", None )
    if not potential_values_per_frame or not polygons_per_frame:
        raise ValueError ( "recolor_mep_surface: this object has no cached MEP data "
                            "(mep_potential_values/mep_polygons) -- this only works for surfaces "
                            "generated as MEP." )

    vmin      = vmin      if vmin      is not None else getattr ( vobject, "mep_vmin", None )
    vmax      = vmax      if vmax      is not None else getattr ( vobject, "mep_vmax", None )
    cmap_name = cmap_name if cmap_name is not None else getattr ( vobject, "mep_cmap_name", "coolwarm" )

    for frame, frame_data in enumerate ( vobject.surface_trajectory ):
        if frame_data is None or "obital_plus" not in frame_data:
            continue
        potential_values = potential_values_per_frame[frame] if frame < len ( potential_values_per_frame ) else None
        polygons_np       = polygons_per_frame[frame]         if frame < len ( polygons_per_frame )         else None
        if potential_values is None or polygons_np is None:
            continue

        vertex_colors = mep_colormap ( potential_values, vmin = vmin, vmax = vmax, cmap_name = cmap_name )
        tri_colors    = vertex_colors[polygons_np].reshape ( -1 ).astype ( np.float32 )

        vertices, _old_colors, indexes, normals = frame_data["obital_plus"]
        frame_data["obital_plus"] = [ vertices, tri_colors, indexes, normals ]

    vobject.mep_vmin      = vmin
    vobject.mep_vmax      = vmax
    vobject.mep_cmap_name = cmap_name


def _pdynamo_array1d_to_numpy ( arr, dtype ):
    """ [EN] 1-D analogue of _pdynamo_array_to_numpy() above, for things
    like potentialProperty.gridValues (a RealArray1D, not an Array2D --
    no .rows/[i,:] row-slicing, just a flat sequence of n scalars). """
    try:
        result = np.asarray ( arr, dtype = dtype )
        if result.ndim == 1 and result.shape[0] == len ( arr ):
            return result
    except Exception:
        pass
    n = len ( arr )
    out = np.empty ( n, dtype = dtype )
    for i in range ( n ):
        out[i] = arr[i]
    return out


def _reconstruct_regular_grid ( pts, vals ):
    """ [EN] THE key optimisation of changelog item 13. Reconstructs a
    proper (values_3d, origin, spacing) regular-grid representation from
    flat, PARALLEL (point, value) pairs whose internal storage order is
    unknown/unverified (pDynamo3's own potentialProperty.gridPoints /
    .gridValues) -- WITHOUT trusting any assumption about that order.

    How: the points are known, by construction (they come from
    generator.DefineGrid()), to lie exactly on a regular 3D grid. So
    instead of guessing how pDynamo3 laid them out internally, we
    discover the grid's own axis structure empirically, straight from the
    coordinates themselves: collect the unique x/y/z coordinate values
    (a true regular grid has exactly nx/ny/nz distinct values per axis,
    evenly spaced), derive origin + spacing from those, then compute each
    point's (i,j,k) grid index directly from its own coordinates and
    scatter its value into the right cell of a fresh (nx,ny,nz) array.
    The RESULT is independent of the order `pts`/`vals` were given in --
    reconstructing from a deliberately SHUFFLED point order was used to
    verify this during development, and reproduced the exact expected
    origin/spacing.

    Once values_3d/origin/spacing exist, _trilinear_interpolate() (the
    same, already-verified function used for external .cube files) can be
    reused directly -- both faster AND more accurate than the nearest-
    neighbour lookup this function replaces (benchmarked: ~0.03s total
    for reconstruction + 30k interpolated queries against a realistic
    ~157k-point grid, versus tens of seconds -- or, with the naive
    chunk_size that shipped briefly, a MemoryError trying to allocate
    14+ GB -- for brute-force nearest-neighbour on the same input).

    Returns None if `pts` does NOT form a complete regular grid (e.g. some
    QC codes only evaluate properties within a cutoff distance of the
    molecule, producing a "hollowed out" or irregular point cloud) -- the
    caller must fall back to _nearest_neighbor_lookup() in that case. """
    xs_u = np.unique ( np.round ( pts[:,0], 6 ) )
    ys_u = np.unique ( np.round ( pts[:,1], 6 ) )
    zs_u = np.unique ( np.round ( pts[:,2], 6 ) )
    nx, ny, nz = len ( xs_u ), len ( ys_u ), len ( zs_u )
    if nx * ny * nz != len ( pts ):
        return None   # not a complete regular box -- caller uses the fallback

    ox, oy, oz = xs_u[0], ys_u[0], zs_u[0]
    dx = ( xs_u[-1] - xs_u[0] ) / ( nx - 1 ) if nx > 1 else 1.0
    dy = ( ys_u[-1] - ys_u[0] ) / ( ny - 1 ) if ny > 1 else 1.0
    dz = ( zs_u[-1] - zs_u[0] ) / ( nz - 1 ) if nz > 1 else 1.0

    ix = np.round ( ( pts[:,0] - ox ) / dx ).astype ( np.int64 )
    iy = np.round ( ( pts[:,1] - oy ) / dy ).astype ( np.int64 )
    iz = np.round ( ( pts[:,2] - oz ) / dz ).astype ( np.int64 )

    values_3d = np.empty ( (nx, ny, nz), dtype = np.float64 )
    values_3d[ix, iy, iz] = vals   # scatter vetorizado -- uma atribuicao so
    return values_3d, (ox, oy, oz), (dx, dy, dz)


def _nearest_neighbor_lookup ( pts, vals, query_points, max_chunk_bytes = 64 * 1024 * 1024 ):
    """ [EN] FALLBACK ONLY, used by build_potential_interpolator() when
    _reconstruct_regular_grid() returns None (points don't form a
    complete regular box). Pure-numpy replacement for
    scipy.spatial.cKDTree -- no assumption about point storage order,
    like the old version, but rewritten to actually be memory-safe:
    the previous version (chunk_size=4000, one (chunk,n,3) tensor per
    batch) tried to allocate 14+ GB and crashed with MemoryError against
    a realistic ~157k-point grid the very first time this ran against
    real data. Fixed by (a) expanding the squared-distance formula
    |a-b|^2 = |a|^2 + |b|^2 - 2*a.b to work with a (chunk, n) MATRIX
    instead of a (chunk, n, 3) TENSOR (an automatic 3x memory cut), and
    (b) choosing chunk_size dynamically so each batch stays under
    max_chunk_bytes (default 64 MB) regardless of how large the grid is,
    instead of a fixed chunk_size that could silently blow up on a larger
    grid than whatever it happened to be tuned against. """
    pts = np.asarray ( pts, dtype = np.float64 )
    vals = np.asarray ( vals, dtype = np.float64 )
    query_points = np.asarray ( query_points, dtype = np.float64 )
    n = pts.shape[0]
    m = query_points.shape[0]
    pts_sq = np.einsum ( 'ij,ij->i', pts, pts )   # |b|^2, (n,) -- calculado uma vez so

    bytes_per_query_row = max ( 1, n * 8 )   # float64
    chunk_size = max ( 1, int ( max_chunk_bytes // bytes_per_query_row ) )

    out = np.empty ( m, dtype = np.float64 )
    for start in range ( 0, m, chunk_size ):
        end = min ( start + chunk_size, m )
        chunk = query_points[start:end]
        chunk_sq = np.einsum ( 'ij,ij->i', chunk, chunk )   # |a|^2, (c,)
        d2 = chunk_sq[:, None] + pts_sq[None, :] - 2.0 * ( chunk @ pts.T )   # (c, n) -- 2D matrix, not 3D tensor
        idx = np.argmin ( d2, axis = 1 )
        out[start:end] = vals[idx]
    return out


def build_potential_interpolator ( potentialProperty ):
    # [EN] Tries _reconstruct_regular_grid() + _trilinear_interpolate()
    # first (fast AND accurate AND makes no assumption about pDynamo3's
    # internal array storage order -- see _reconstruct_regular_grid()'s
    # docstring for the full reasoning and benchmark numbers). Falls back
    # to _nearest_neighbor_lookup() only if the potential grid turns out
    # not to be a complete regular box (rare, but not impossible -- some
    # QC codes only evaluate grid properties within a cutoff of the
    # molecule).
    """ Recebe o QCGridProperty bruto do potencial (generator.GetProperty(tag),
    ANTES de virar isosuperficie -- precisa ter .gridPoints/.gridValues) e
    devolve uma funcao que avalia o potencial em qualquer ponto 3D. """
    pts_np  = _pdynamo_array_to_numpy ( potentialProperty.gridPoints, np.float64 )
    vals_np = _pdynamo_array1d_to_numpy ( potentialProperty.gridValues, np.float64 )

    reconstructed = _reconstruct_regular_grid ( pts_np, vals_np )
    if reconstructed is not None:
        values_3d, origin, spacing = reconstructed

        def evaluate ( query_points ):
            """ query_points: array (m,3), nas MESMAS unidades do grid (Bohr,
            que e a unidade nativa do pDynamo -- ver surface_parser, que so
            converte pra Angstrom na hora de montar o buffer de exibicao). """
            return _trilinear_interpolate ( values_3d, origin, spacing, query_points )
    else:
        dprint ( "build_potential_interpolator: the potential grid does not form a "
                "complete regular box -- using nearest-neighbor (slower) "
                "instead of trilinear interpolation." )

        def evaluate ( query_points ):
            return _nearest_neighbor_lookup ( pts_np, vals_np, query_points )

    return evaluate


def surface_parser_mep ( surface, vertex_colors ):
    # [EN] VECTORISED (changelog item 13), same technique as surface_parser()
    # above -- see that function's comment for the benchmark numbers.
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
    if skipped:
        dprint ( "surface_parser_mep: {} triangulo(s) fantasma (aresta anormal) descartado(s)".format ( skipped ) )

    polygons_np      = _pdynamo_array_to_numpy ( polygons, np.int64 )[valid_mask]
    vertices_np      = _pdynamo_array_to_numpy ( vertices, np.float64 )
    normals_np       = _pdynamo_array_to_numpy ( smooth_normals, np.float64 )
    vertex_colors_np = np.asarray ( vertex_colors, dtype = np.float64 )

    tri_verts  = vertices_np[polygons_np] / 1.889725989   # (n_valid_tri, 3, 3) -- Bohr -> Angstrom
    tri_norms  = normals_np[polygons_np]                  # (n_valid_tri, 3, 3)
    tri_colors = vertex_colors_np[polygons_np]            # (n_valid_tri, 3, 3)

    vertices_out = tri_verts.reshape ( -1 ).astype ( np.float32 )
    normals_out  = tri_norms.reshape ( -1 ).astype ( np.float32 )
    colors_out   = tri_colors.reshape ( -1 ).astype ( np.float32 )

    # one index per vertex (not per component/float) -- unlike the
    # surface_parser original, que gera indexes com 3x mais entradas do
    # que vertices de verdade existem no buffer (ver nota no README).
    indexes_out = np.arange ( vertices_out.shape[0] // 3, dtype = np.uint32 )
    # [EN] Also returns polygons_np (the mask-filtered triangle index
    # array used just above, tri_colors = vertex_colors_np[polygons_np])
    # -- needed by the caller to cache it alongside the raw, pre-colormap
    # potential_values, so a LATER change to vmin/vmax/colormap (see
    # _surf_setup in treeview_menu.py) can redo just this cheap
    # expansion step (mep_colormap() + vertex_colors_np[polygons_np]) --
    # without re-running MakeVertexNormalsFromPolygonalNormals(),
    # _compute_valid_polygon_mask(), or anything that touches pDynamo at
    # all. Was a 4-tuple return before; every call site needed updating.
    return vertices_out, colors_out, indexes_out, normals_out, polygons_np


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
        # branch 'mep' -- except that the two cubes (density and potential)
        # vem de arquivos externos em vez do grid do pDynamo.
        potential_cube     = read_cube_file ( potential_path )
        evaluate_potential = build_potential_interpolator_from_cube ( potential_cube )
        # [EN] vectorised conversion (changelog item 13) -- was a
        # per-vertex, per-component Python loop before.
        verts_bohr = _pdynamo_array_to_numpy ( surface.vertices, np.float64 )
        potential_values = evaluate_potential ( verts_bohr )
        vertex_colors    = mep_colormap ( potential_values, vmin = mep_vmin, vmax = mep_vmax,
                                           cmap_name = mep_cmap_name )
        vertices, colors, indexes, normals, polygons_np = surface_parser_mep ( surface, vertex_colors )
        # [EN] cached alongside the mesh so vmin/vmax/colormap can be
        # changed later (see _surf_setup in treeview_menu.py) without
        # re-reading the .cube files or re-running marching cubes --
        # see the matching cache entries in the 'mep' branch of
        # generate_grid_parallel() below for the full reasoning.
        return { 'obital_plus'              : [ vertices, colors, indexes, normals ],
                 'mep_raw_potential_values'  : potential_values,
                 'mep_raw_polygons'          : polygons_np }
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
    _t_start = time.perf_counter ( )   # TEMPORARY DEBUG -- see print at the end of the function

    # [EN] Early-return guard, added for the "External" cube-import surface
    # type (changelog item 9). Every other branch below unconditionally
    # calls apply_coords_to_system() / system.Energy() /
    # QCGridPropertyGenerator.FromSystem(system) a few lines further down,
    # all of which require a real, live pDynamo QC system -- but an
    # externally-supplied .cube file has no such system (system/coords are
    # simply None in that job tuple). Must return BEFORE that setup code
    # runs, not after.
    if _type == 'external_cube':
        # external cube (.cube, e.g. ORCA via orca_plot) -- does not need
        # system/coords/QCGridPropertyGenerator nenhum, so leitura de
        # file. Branches BEFORE the QC system setup just below, which
        # exigiria um "system" de verdade (None aqui, ja que essa entrada
        # does not come from a pDynamo QC calculation).
        return _generate_external_cube_surface ( parameters )

    _GridSpacing   = parameters['_GridSpacing']
    _OrbitalTag    = parameters['_OrbitalTag']
    _isovalue      = parameters['_isovalue']  
    _IsosurfaceTag = parameters['_IsosurfaceTag']
    _mep_vmin      = parameters.get ( 'mep_vmin' )   # None = automatico (percentil)
    _mep_vmax      = parameters.get ( 'mep_vmax' )
    _mep_cmap_name = parameters.get ( 'mep_cmap_name', 'coolwarm' )
    _mep_pot_spacing = parameters.get ( 'mep_pot_spacing' )   # None = automatico (2.5x _GridSpacing)
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
    _t_grid0 = time.perf_counter ( )   # DEBUG TEMPORARIO
    generator.DefineGrid    ( gridSpacing = _GridSpacing ) # . Some value in atomic units - e.g. 0.2
    dprint ( "DEBUG TIMING: DefineGrid (gridSpacing={}) levou {:.3f} s".format ( _GridSpacing, time.perf_counter() - _t_grid0 ) )
    

    
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
        # 1. Mesh geometry from the DENSITY isosurface
        #    (the isovalue in the entry_isovalue field now means
        #    "density isovalue" in this mode -- ~0.002-0.02 a.u. is usually
        #    aproximar bem o contorno de van der Waals).
        _t = time.perf_counter ( )   # DEBUG TEMPORARIO
        generator.GridDensity ( tag = 'density_mep' )
        dprint ( "DEBUG TIMING: GridDensity levou {:.3f} s".format ( time.perf_counter() - _t ) ); _t = time.perf_counter ( )
        generator.Isosurface  ( 'density_mep', _isovalue, tag = _IsosurfaceTag )
        dprint ( "DEBUG TIMING: Isosurface (density) levou {:.3f} s".format ( time.perf_counter() - _t ) ); _t = time.perf_counter ( )
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        density_iso = surfaceProperty.isosurface
        dprint ( "DEBUG TIMING: n_vertices={} n_triangulos={}".format ( density_iso.vertices.rows, density_iso.polygons.rows ) )

        # 2. Grid de POTENCIAL bruto (os valores a mapear na malha acima).
        #    Tag propria, distinta de 'density_mep' e de _IsosurfaceTag --
        #    ver o bug de colisao de tags que corrigimos no branch 'potential'.
        # [EN] KEY OPTIMISATION (changelog item 14): GridPotential() was
        # measured to be ~180x slower than GridDensity() on the SAME grid
        # (confirmed inherent to the method -- see QCModelBase.py's
        # GridPointPotentials, an O(n_basis^2)-per-point calculation, vs
        # GridPointDensities' O(n_basis)-per-point). Since MEP already
        # interpolates the potential onto the density mesh via real
        # trilinear interpolation (_reconstruct_regular_grid +
        # _trilinear_interpolate), the potential grid does NOT need to
        # share the density grid's (possibly very fine) spacing -- a much
        # coarser grid, just for this step, is still interpolated
        # smoothly onto the finer density-surface vertices. Confirmed
        # safe to call DefineGrid() again here, after GridDensity(): each
        # QCGridProperty stores its OWN grid/gridPoints reference at the
        # moment it is computed (see QCGridProperties.py), so the density
        # isosurface computed above keeps referring to the FINE grid it
        # was actually computed on, unaffected by redefining self.grid
        # for the potential step below.
        _potential_spacing = _mep_pot_spacing if _mep_pot_spacing is not None else _GridSpacing * 2.5
        generator.DefineGrid ( gridSpacing = _potential_spacing )
        generator.GridPotential ( tag = 'potential_mep' )
        dprint ( "DEBUG TIMING: GridPotential (spacing={:.3f}, {}x coarser than density) took {:.3f} s".format (
                _potential_spacing, _potential_spacing / _GridSpacing, time.perf_counter() - _t ) ); _t = time.perf_counter ( )
        potentialProperty  = generator.GetProperty ( 'potential_mep' )
        evaluate_potential = build_potential_interpolator ( potentialProperty )
        dprint ( "DEBUG TIMING: build_potential_interpolator (inclui reconstrucao do grid) levou {:.3f} s".format ( time.perf_counter() - _t ) ); _t = time.perf_counter ( )

        # 3. Avalia o potencial em cada vertice ORIGINAL da malha de densidade,
        #    em Bohr (mesma unidade do grid do pDynamo -- a conversao pra
        #    Angstrom so acontece dentro de surface_parser_mep).
        # [EN] vectorised conversion (changelog item 13) -- was a
        # per-vertex, per-component Python loop before.
        verts_bohr = _pdynamo_array_to_numpy ( density_iso.vertices, np.float64 )

        potential_values = evaluate_potential ( verts_bohr )
        dprint ( "DEBUG TIMING: evaluate_potential (interpolacao) levou {:.3f} s".format ( time.perf_counter() - _t ) ); _t = time.perf_counter ( )
        vertex_colors    = mep_colormap ( potential_values, vmin = _mep_vmin, vmax = _mep_vmax,
                                           cmap_name = _mep_cmap_name )
        dprint ( "DEBUG TIMING: mep_colormap levou {:.3f} s".format ( time.perf_counter() - _t ) ); _t = time.perf_counter ( )

        vertices, colors, indexes, normals, polygons_np = surface_parser_mep ( density_iso, vertex_colors )
        dprint ( "DEBUG TIMING: surface_parser_mep levou {:.3f} s".format ( time.perf_counter() - _t ) )
        orbital_iso['obital_plus'] = [vertices, colors, indexes, normals]
        # [EN] Cached here (alongside the mesh itself) so a LATER change
        # to the MEP colour scale (vmin/vmax/colormap -- see _surf_setup
        # in treeview_menu.py, and the user's original request: "quero
        # mexer no vmin e vmax do MEP... todo o calculo tem que ser
        # refeito, isso nao e bom") only needs to redo the cheap part:
        # mep_colormap(potential_values, ...) + vertex_colors[polygons_np]
        # -- NOT GridDensity(), Isosurface(), GridPotential(),
        # build_potential_interpolator(), or MakeVertexNormalsFromPolygonalNormals(),
        # all of which are the genuinely expensive steps above (see the
        # DEBUG TIMING prints throughout this branch -- GridPotential
        # alone was measured at ~180x the cost of GridDensity on the
        # same grid). potential_values is per ORIGINAL density-surface
        # vertex (pre-colormap, pre-triangle-expansion); polygons_np is
        # the mask-filtered triangle index array surface_parser_mep()
        # used to expand per-vertex colours into per-triangle-corner
        # colours (vertex_colors_np[polygons_np]) -- both are exactly
        # what's needed to redo just that expansion with a new colour
        # scale, cheaply, in the main process, without touching pDynamo
        # or the multiprocessing pool at all.
        orbital_iso['mep_raw_potential_values'] = potential_values
        orbital_iso['mep_raw_polygons']         = polygons_np
    
    
    
    
    
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
    
    dprint ( "DEBUG TIMING: generate_grid_parallel TOTAL levou {:.3f} s (tipo: {})".format ( time.perf_counter() - _t_start, _type ) )
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

    # [EN] Captures the LCAO (linear combination of atomic orbitals)
    # coefficient matrix for this frame -- user's request: "verifique se
    # o pDynamo permite capturar qual a combinacao linear que gerou o
    # orbital molecular" (confirmed: yes -- orbitalsP.orbitals is exactly
    # this, a (numberBasisFunctions, numberOrbitals) matrix; see
    # pMolecule/QCModel/QCOrbitals.py in pDynamo3's own source, where
    # MakeFromFock() fills it in directly from EigenPairs() on the Fock
    # matrix -- each COLUMN j is orbital j's coefficient vector over the
    # atomic-orbital basis). Converted to numpy immediately, in THIS
    # worker process, rather than relying on the returned `system` object
    # still having a fully-usable orbitalsP.orbitals after being pickled
    # back across the multiprocessing boundary.
    #
    # Also captures, so individual coefficients can be attributed back to
    # a specific QC atom later: system.qcState.orbitalBases.
    # centerFunctionPointers (basis-function index boundaries per QC
    # atom -- the same array pDynamo3's own Mulliken population-analysis
    # code uses for exactly this purpose, see
    # MullikenMultipoleEvaluator.pyx) and system.qcState.atomicNumbers
    # (converted to element symbols via pScientific.PeriodicTable, in the
    # SAME order as centerFunctionPointers).
    try:
        orbitals_matrix = _pdynamo_array_to_numpy ( orbitalsP.orbitals, np.float64 )
        center_function_pointers = _pdynamo_array1d_to_numpy (
            system.qcState.orbitalBases.centerFunctionPointers, np.int64 )
        atomic_numbers = _pdynamo_array1d_to_numpy ( system.qcState.atomicNumbers, np.int64 )
        from pScientific import PeriodicTable
        atom_symbols = [ PeriodicTable.Symbol ( int ( z ) ) for z in atomic_numbers ]
    except Exception as e:
        # Nao deixa a falta/erro da parte de LCAO quebrar o calculo de
        # energias/orbitais em si, que ja funcionava antes disso existir.
        dprint ( "generate_wavefunction_parallel: could not capture LCAO data ({}).".format ( e ) )
        orbitals_matrix = None
        center_function_pointers = None
        atom_symbols = None

    return orbitals, system, generator, orbitals_matrix, center_function_pointers, atom_symbols
    
    
    
    
    
    '''
    for i in range(len(orbitals)):
        reverse_index = -i-1 #- len(orbitals)
        dprint(reverse_index, orbitals[reverse_index ])
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


