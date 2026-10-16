#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  molecular_surface.py
#
#  EDTSurf-equivalent molecular surface generator (Van der Waals / Solvent-
#  Accessible / Solvent-Excluded surfaces) -- user's own explicit request:
#  "a representacao de superficie, tal como temos em
#  https://zhanggroup.org/EDTSurf/". Pure numpy + pDynamo3's own compiled
#  marching-cubes routine, NO scipy dependency (see "no scipy" note below),
#  no GTK dependency either (mirrors util/colormaps.py's own importable-
#  anywhere convention) -- the GL/GTK side (SurfaceRepresentation, the
#  treeview menu's own pre-generation setup dialog) is wired up separately
#  and just calls build_surface_mesh() below.
#
#  [EN] Algorithm summary:
#    - VWS (Van der Waals surface): the isosurface at 0 of the ANALYTIC
#      signed field F(p) = min_i(|p - atom_i| - vdw_radius_i) -- a plain
#      union-of-spheres, no distance-transform needed.
#    - SAS (Solvent-Accessible Surface): the SAME field, but with each
#      atom's radius inflated by `probe_radius` first (the classic "roll a
#      probe sphere, trace its CENTER" surface).
#    - SES (Solvent-Excluded / Connolly surface): the one type with no
#      simple analytic per-atom formula -- EDTSurf's own actual algorithmic
#      contribution. Built via the classic two-pass distance-transform
#      construction: a point is truly solvent-EXCLUDED (part of the SES
#      solid) iff it is SAS-occupied AND every point within `probe_radius`
#      of it is ALSO SAS-occupied (i.e. no probe-sphere placement can ever
#      reach it) -- everything else (deep pockets a probe can't reach,
#      inter-atom cavities too small for the probe) gets correctly carved
#      back OUT of the SAS solid. See _bounded_edt_from_true_to_false()
#      below for the (scipy-free) implementation.
#
#  [EN] NO SCIPY (user's own explicit request, 2026-09-23: "tente eliminar
#  a dependencia do scipy e use apenas numpy"). Two scipy calls this module
#  used to make, and what replaced them:
#    - scipy.spatial.cKDTree (whole-grid nearest-atom field + nearest-atom
#      vertex colouring) -> replaced by a LOCAL SPLAT (each atom only ever
#      updates its own small sub-block of the grid, an elementwise
#      np.minimum against whatever's already there -- see
#      _signed_union_of_spheres_field()) for the field, and a plain
#      CHUNKED BRUTE-FORCE nearest-atom search (see _nearest_atom_colors())
#      for vertex colouring -- exact, not approximate (the earlier
#      cKDTree(k=12) version was itself only an approximation), and cheap
#      because pDynamo3's own marching-cubes vertex array is ALREADY a
#      welded/shared-vertex mesh (confirmed empirically: a synthetic sphere
#      test showed every single output vertex shared by >1 triangle, ~6 on
#      average) -- so colouring only needs to run once per UNIQUE vertex,
#      not once per triangle corner.
#    - scipy.ndimage.distance_transform_edt (SES carving) -> replaced by
#      _bounded_edt_from_true_to_false(), an exact (within a small
#      truncation radius) Euclidean distance transform computed via a
#      brute-force scan over every integer grid offset within that radius
#      (a handful of voxels at the spacings this module uses -- e.g. ~30
#      offsets for probe_radius=1.4/spacing=0.8), each one a single
#      vectorised numpy comparison over the whole grid. This is exact
#      (unlike a chamfer/approximate distance transform) precisely because
#      SES only ever needs distances up to `probe_radius` (the boolean
#      "is this point solvent-excluded" decision) or up to a couple of
#      voxels (the smooth field marching cubes interpolates near the
#      zero-crossing) -- never a TRUE whole-grid distance transform, so
#      truncating is free correctness-wise, not an approximation.
#
#  MarchingCubes_Isosurface3D (from pDynamo3's pScientific.Surfaces) is
#  fully generic -- it only needs a RegularGrid + a 3D scalar array + an
#  isovalue, no QC/wavefunction dependency (confirmed by this app's own
#  existing cube_to_pdynamo_surface() in gui/windows/analysis/
#  surface_analysis_window.py, which builds a grid from an arbitrary
#  external .cube file the same way). This module reuses that exact same
#  engine, plus two already-debugged helpers from that same file
#  (_pdynamo_array_to_numpy, _compute_valid_polygon_mask -- see below)
#  rather than re-deriving them.
#
#  [EN] REAL BUG, FOUND TWICE, FIXED PROPERLY THE SECOND TIME (2026-09-23):
#  pDynamo3's OWN compiled marching-cubes C routine (MarchingCubes.c)
#  silently CORRUPTS its output once a mesh needs more vertices than its
#  own initial buffer allocation (`_VertexCount0` = 10000, but the REAL
#  initial size is `min(3*nCubes, _VertexCount0)` -- i.e. it can trigger
#  at vertex counts well below 10000 too, depending on how many grid cubes
#  the isosurface crosses) -- confirmed by reading that C source directly.
#  This is the EXACT SAME bug already found, root-caused, and fixed for
#  this app's own QC (orbital/density/potential) surface tool -- see
#  _compute_valid_polygon_mask() in surface_analysis_window.py, whose own
#  docstring has the full derivation (corrupted vertices land at the
#  REGULAR GRID's own "midPointLower" origin corner, not the world
#  origin -- easy to misdiagnose as "connecting to the origin" when the
#  molecule itself isn't centred near world (0,0,0)).
#
#  This module's FIRST attempt at a fix (now removed) checked for an EXACT
#  match against that one analytically-predicted corrupted-vertex position
#  -- which turned out to be insufficient: the VWS surface type kept
#  showing the same "points connecting to one corner" artifact the user
#  reported, because the real corruption trigger depends on `nCubes`
#  (different per surface type/grid shape), not a single universal vertex
#  count, so an exact-point match can miss cases the QC tool's own filter
#  already handled. Fixed properly by reusing THAT already-proven filter
#  instead: _compute_valid_polygon_mask()'s GEOMETRIC (not positional)
#  approach -- discard any triangle whose longest edge is a wild outlier
#  relative to the mesh's own median edge length. General by construction,
#  independent of grid shape/vertex-count thresholds, and already battle-
#  tested on this exact class of bug.
#
#  [EN] SMOOTH NORMALS (point II of the user's 2026-09-23 follow-up
#  request, "os poligonos estao com faces bem definidas, seria
#  interessante combinar as normais"): investigated and found this was
#  NEVER actually a data problem -- surface.MakeVertexNormalsFromPolygonalNormals()
#  (called below) already produces correctly-AVERAGED per-vertex normals,
#  confirmed empirically (a synthetic sphere: 1854 unique vertices, every
#  one shared by >1 triangle, each normal genuinely being the average of
#  its ~6 adjacent face normals -- pDynamo3's own vertices/polygons arrays
#  ARE a real welded/shared-vertex mesh, not per-triangle-duplicated). The
#  faceted look was caused purely by SurfaceRepresentation's own
#  `smooth_shading` flag defaulting to False (flat) -- the fix belongs
#  entirely in the GL/GTK layer (treeview_menu.py calls
#  representation.set_shading_mode("smooth") right after creating a
#  molecular-surface representation), not in this module.
#
#  [EN] Unit note: unlike surface_analysis_window.py's own surface_parser()
#  (which divides vertices by 1.889725989, Bohr -> Angstrom, because QC
#  grids are natively in Bohr), THIS module's grid is built directly from
#  vismol atom coordinates, which are already in Angstrom -- no conversion
#  is applied anywhere here. Do not reuse surface_parser() itself for this
#  module's output; it would wrongly shrink the mesh ~1.89x.
# ============================================================================
import numpy as np

from pScientific.Geometry3 import RegularGrid
from pScientific.Arrays    import Array, Reshape, RealArrayND
from pScientific.Surfaces  import MarchingCubes_Isosurface3D
from pScientific.Surfaces.PolygonalSurface import PolygonalSurface


DEFAULT_PROBE_RADIUS = 1.4    # Angstrom -- standard water probe (EDTSurf's own default)
DEFAULT_GRID_SPACING = 0.8    # Angstrom per voxel -- see the pDynamo3 vertex-buffer bug
                               # note above: coarse enough to stay clear of it for most
                               # small/medium proteins; the retry loop in build_surface_mesh()
                               # is the real safety net for anything bigger.
_LOCAL_SPLAT_MARGIN_CELLS = 2   # extra empty voxels of padding around each atom's own
                                 # local field sub-block (see _signed_union_of_spheres_field)
_BOUNDING_BOX_MARGIN_CELLS = 3   # extra empty voxels of padding around the WHOLE molecule
                                  # (so the isosurface never gets clipped right at the box edge)

# [EN] pDynamo3 marching-cubes corruption workaround (see top-of-file note):
# if more than this fraction of output triangles get discarded by the
# (proven, geometric) ghost-triangle filter, treat the WHOLE result as
# unusable and retry coarser -- a couple of discarded triangles is normal/
# expected numerical noise, this threshold only catches genuine corruption.
_MAX_GHOST_TRIANGLE_FRACTION = 0.01
_MAX_RESOLUTION_RETRIES = 5
_GRID_SPACING_RETRY_GROWTH = 1.35

_COLOR_CHUNK_SIZE = 4096   # vertices per chunk for the brute-force nearest-atom colour search


def _isosurface_from_grid ( grid, dataND, isovalue ):
    """ [EN] PERFORMANCE FIX (profiling found this was ~35% of this
    module's total per-frame time after the other fixes above --
    confirmed via cProfile on a real crambin.pdb surface). Calls
    pDynamo3's own `MarchingCubes_Isosurface3D()`, but temporarily no-ops
    `PolygonalSurface.MakePolygonNormals` for the duration of that ONE
    call. Confirmed by reading MarchingCubes.pyx directly: the wrapper
    ALREADY does `surface = PolygonalSurface.FromArrays(...) ; surface.
    MakePolygonNormals() ; return surface` -- and the C step underneath
    (CMarchingCubes_Isosurface3D) already computes `polygonNormals`
    itself as one of its own output arrays, so that Python-level call is
    PURELY REDUNDANT unless the caller actually reads `surface.
    polygonNormals`/`surface.vertexNormals` afterward -- this module
    never does (see _vertex_normals_from_mesh() above, our own much
    faster numpy replacement), so the ~0.2s it costs (a per-triangle,
    pure-Python loop of Vector3 method calls, see PolygonalSurface.
    MakePolygonNormals()'s own source) was pure waste for this module's
    own use case.

    Safe: this only patches the CLASS method for the exact duration of
    this synchronous call (try/finally, restored immediately after,
    before this function returns) -- no other code runs in between in
    this single-threaded codebase, and each of this module's own
    multiprocessing.Pool WORKER processes (see build_surface_trajectory()
    below) has its own independent copy of the imported module, so this
    can never affect a concurrent call elsewhere (e.g. the QC surface
    tool, which DOES need real polygon/vertex normals and is never
    patched by this). Does not touch any file on disk -- a pure runtime
    patch, gone the instant this function returns. """
    original_make_polygon_normals = PolygonalSurface.MakePolygonNormals
    PolygonalSurface.MakePolygonNormals = lambda self: None
    try:
        return MarchingCubes_Isosurface3D ( grid, dataND, isovalue )
    finally:
        PolygonalSurface.MakePolygonNormals = original_make_polygon_normals


COLOR_MODES = ( "atom", "chain", "carbon" )   # see _resolve_atom_colors()'s own docstring

# [EN] Fixed cycle used by "chain" coloring -- distinct, easy-to-tell-apart
# hues, deliberately NOT reusing util/colorpalette.py's own CPK-style
# per-ELEMENT palette (that one colors atoms by what they ARE, chemically;
# this one just needs N visually distinct buckets for however many chains
# a structure happens to have, so a short hand-picked cycle is simpler and
# more distinguishable than trying to repurpose an element palette).
_CHAIN_COLOR_CYCLE = [
    [ 0.90, 0.30, 0.30 ], [ 0.30, 0.55, 0.90 ], [ 0.35, 0.75, 0.35 ], [ 0.95, 0.70, 0.20 ],
    [ 0.65, 0.35, 0.85 ], [ 0.30, 0.80, 0.80 ], [ 0.90, 0.45, 0.70 ], [ 0.55, 0.55, 0.55 ],
]


def _resolve_atom_colors ( vismol_object, atom_ids, atoms, color_mode ):
    """ [EN] Returns a (n_atoms, 3) float32 RGB array, ONE ROW PER ATOM
    (same order as `atom_ids`/`atoms`), for whichever `color_mode` was
    requested -- resolved ONCE here, on the caller's process, BEFORE
    atoms cross into the plain-numpy/picklable per-frame pipeline (see
    _build_surface_mesh_from_atoms()/_generate_surface_frame_job() above)
    -- so nothing downstream (nearest-atom sampling in
    _nearest_atom_colors()) ever needs to know a "color mode" concept
    exists at all; it always just samples a plain per-atom RGB array,
    the exact same way regardless of mode. This is also what keeps the
    door open for the user's own explicitly requested future mode
    ("eletrostatico" -- electrostatic potential, not implemented yet,
    needs per-atom partial charges this app doesn't universally assign):
    it would be one more branch here, producing one more (n_atoms,3)
    array (e.g. via a diverging colormap over per-atom charge), with
    ZERO changes needed anywhere else in this module.

    User's own explicit request, 2026-09-23: "by chain, usar as cores
    dos carbonos do sistema, e outros (deixe em aberto para receber o
    potencial eletrostatico depois)".

    "atom" (default, previous/only behaviour): each atom keeps its own
    existing EasyHybrid display color (`vismol_object.colors`, already
    populated for every other representation -- CPK by element, unless
    overridden).
    "chain": one distinct color per chain (`atom.chain.name`, `_CHAIN_
    COLOR_CYCLE` above, cycling if there are more chains than colors);
    atoms with no chain (`atom.chain is None`) share a single "" bucket.
    "carbon": per-atom, same as the standard CPK palette for every
    element EXCEPT carbon, which uses this object's own "Reference
    Color" (`vismol_object.color_palette['C']`, the same per-system
    carbon-color override already used app-wide, see e.g. treeview_
    menu.py's own 'Reference Color' menu item / change_reference_color()
    in main_window.py, to tell different loaded molecules apart) --
    NOT a single uniform color for the whole surface (the user's own
    clarification, 2026-09-23: "os outros atomos (nao carbonos) tenham
    as cores tal como sao na paleta"). Reads straight from `vismol_
    object.color_palette` (the {symbol: RGB} dict itself, already
    carrying this exact carbon override -- see PeriodicTable.
    get_color_palette()) rather than `vismol_object.colors` (point
    "atom" above), since the latter is whatever this object's CURRENT
    on-screen coloring happens to be (could include other, unrelated
    per-atom overrides) -- this mode is meant to be a clean, deterministic
    "CPK + this system's own carbon color", regardless of that.

    [EN] REAL BUG FOUND AND FIXED (2026-09-24, right after the "not a
    single uniform color" fix above -- the user's own "voce nao
    entendeu" was this, not a logic misunderstanding): confirmed live
    that `main_window.py`'s own `change_reference_color()` -- the ACTUAL
    code behind this app's "Reference Color" menu item, i.e. the normal,
    common way a real user sets "their respective system"'s carbon color
    -- does `system.e_color_palette['C'] = new_color` with WHATEVER it's
    given, unnormalized (a `Gdk.RGBA`-derived 4-tuple in the real UI
    flow, confirmed live: `(1.0, 0.0, 1.0, 1.0)`), unlike every OTHER
    entry in that same dict (plain 3-component RGB, from PeriodicTable).
    Building a numpy array that mixes a 4-tuple 'C' entry with 3-tuple
    entries for every other element raised `ValueError: setting an
    array element with a sequence` -- i.e. `color_mode="carbon"` was
    silently BROKEN (a hard crash) for any system whose reference color
    was ever actually changed via that menu item, which is precisely the
    realistic case this whole feature exists for. Fixed by normalizing
    every palette lookup to its first 3 components via `_rgb3()` below,
    rather than trusting palette values to always be clean RGB triples. """
    if color_mode == "atom":
        return np.asarray ( vismol_object.colors, dtype = np.float32 )[ atom_ids ]

    if color_mode == "carbon":
        def _rgb3 ( value, fallback = ( 0.5, 0.5, 0.5 ) ):
            arr = np.asarray ( value, dtype = np.float32 ).reshape ( -1 )
            return arr[ :3 ] if arr.size >= 3 else np.asarray ( fallback, dtype = np.float32 )

        palette = vismol_object.color_palette
        carbon_rgb = _rgb3 ( palette.get ( "C", ( 0.5, 0.5, 0.5 ) ) )
        return np.array (
            [ carbon_rgb if atom.symbol == "C" else _rgb3 ( palette.get ( atom.symbol, carbon_rgb ) ) for atom in atoms ],
            dtype = np.float32,
        )

    if color_mode == "chain":
        chain_names = [ ( atom.chain.name if atom.chain is not None else "" ) or "" for atom in atoms ]
        distinct = sorted ( set ( chain_names ) )
        palette = { name: _CHAIN_COLOR_CYCLE[ i % len ( _CHAIN_COLOR_CYCLE ) ] for i, name in enumerate ( distinct ) }
        return np.array ( [ palette[ name ] for name in chain_names ], dtype = np.float32 )

    raise ValueError ( "unknown color_mode {!r} -- expected one of {}".format ( color_mode, COLOR_MODES ) )


def _gather_atoms ( vismol_object, frame = None, color_mode = "atom" ):
    """ Returns (coords (n,3) float64, vdw_radii (n,) float64, colors (n,3)
    float32) for every atom in `vismol_object`, index-aligned with each
    other -- `vismol_object.atoms` is a plain {0: Atom, 1: Atom, ...} dict
    (confirmed via vismol_object.py's own define_molecules() docstring).
    `colors` comes from _resolve_atom_colors() above -- see its own
    docstring for what `color_mode` ("atom"/"chain"/"carbon") means. """
    atom_ids = sorted ( vismol_object.atoms.keys ( ) )
    atoms = [ vismol_object.atoms[i] for i in atom_ids ]
    coords = np.array ( [ atom.coords ( frame ) for atom in atoms ], dtype = np.float64 )

    periodic_table = vismol_object.vm_session.periodic_table
    radii = np.array (
        [ periodic_table.get_vdw_radius ( symbol = atom.symbol ) for atom in atoms ],
        dtype = np.float64,
    )

    colors = _resolve_atom_colors ( vismol_object, atom_ids, atoms, color_mode )
    return coords, radii, colors


def _grid_bounds ( coords, max_radius, spacing ):
    """ Axis-aligned bounding box (world coords) padded by the largest
    sphere that will ever be splatted (max atom radius, already including
    any probe-radius inflation the caller applied) plus a few empty grid
    cells, so the isosurface never gets clipped right at the box edge. """
    pad = max_radius + _BOUNDING_BOX_MARGIN_CELLS * spacing
    lower = coords.min ( axis = 0 ) - pad
    upper = coords.max ( axis = 0 ) + pad
    dims = np.ceil ( ( upper - lower ) / spacing ).astype ( np.int64 ) + 1
    return lower, ( int ( dims[0] ), int ( dims[1] ), int ( dims[2] ) )


def _signed_union_of_spheres_field ( coords, radii, lower, dims, spacing ):
    """ F(p) = min_i(|p - atom_i| - radii_i) over the whole grid -- the
    ANALYTIC signed distance to a union of spheres (negative inside,
    positive outside, exact zero-crossing at the true sphere surfaces, no
    voxel-resolution artifacts).

    [EN] LOCAL SPLAT, no scipy/KD-tree (see this module's own top-of-file
    "NO SCIPY" note): each atom only ever touches its OWN small local
    sub-block of the grid (its own radius + a couple of empty margin
    cells), updated via a single vectorised np.minimum against whatever
    is already there -- O(atoms * local_voxels), independent of total
    grid size, and exact (not a k-nearest-neighbour approximation like
    the cKDTree version this replaced). Voxels no atom ever reaches keep
    a large constant sentinel value (unambiguously "outside", far from
    any real zero-crossing, so it can never itself produce a spurious
    isosurface triangle -- the marching-cubes "vertex-buffer corruption"
    artifact this module already guards against, see top-of-file note,
    was root-caused to the compiled library's own buffer-growth bug, not
    to field flatness, despite an earlier version of this module briefly
    suspecting the latter). """
    nx, ny, nz = dims
    field = np.full ( ( nx, ny, nz ), 1.0e6, dtype = np.float64 )
    lower = np.asarray ( lower, dtype = np.float64 )
    dims_arr = np.array ( [ nx, ny, nz ], dtype = np.int64 )

    for center, radius in zip ( coords, radii ):
        margin = radius + _LOCAL_SPLAT_MARGIN_CELLS * spacing
        lo = np.maximum ( np.floor ( ( center - margin - lower ) / spacing ).astype ( np.int64 ), 0 )
        hi = np.minimum ( np.ceil  ( ( center + margin - lower ) / spacing ).astype ( np.int64 ) + 1, dims_arr )
        if np.any ( hi <= lo ):
            continue

        # [EN] 2026-09-25 BUG FIX (user's own report: the generated
        # surface was "close, but not exactly" centred on the molecule --
        # confirmed to be a real, systematic half-voxel offset, not a
        # unit-conversion error): pDynamo3's own RegularGrid does NOT map
        # grid index i to physical position `lower + i*binSize` -- it maps
        # it to `midPointLower + i*binSize`, where `midPointLower = lower
        # + 0.5*binSize` (confirmed directly in pDynamo3's own C source,
        # RegularGrid.c/Coordinates3.c: "midPointLower = lower + 0.5e+00 *
        # gridSize", used by every index->coordinate conversion pDynamo3
        # itself performs, including MarchingCubes_Isosurface3D's own
        # vertex output). This module's own top-of-file note ALREADY names
        # `midPointLower` as a real RegularGrid concept (in the context of
        # the unrelated vertex-buffer-corruption bug), but this specific
        # sampling loop was still using the raw, un-shifted `lower` value
        # -- so the analytic sphere-union field sampled here was offset by
        # exactly HALF A GRID CELL, in all 3 axes, from the physical
        # position pDynamo3's own marching cubes later reports/interprets
        # that same sample as living at. `RegularGrid.FromDimensionData()`
        # below is UNCHANGED and still correctly passed the raw `lower` --
        # it computes its OWN midPointLower internally from that; only
        # THIS sampling loop needed the explicit +0.5*spacing shift, so
        # what we sample here and what pDynamo3 later reports agree.
        xs = lower[0] + 0.5 * spacing + np.arange ( lo[0], hi[0] ) * spacing
        ys = lower[1] + 0.5 * spacing + np.arange ( lo[1], hi[1] ) * spacing
        zs = lower[2] + 0.5 * spacing + np.arange ( lo[2], hi[2] ) * spacing
        gx, gy, gz = np.meshgrid ( xs, ys, zs, indexing = "ij" )
        local_field = np.sqrt ( ( gx - center[0] ) ** 2 + ( gy - center[1] ) ** 2 + ( gz - center[2] ) ** 2 ) - radius

        sub = field[ lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2] ]
        np.minimum ( sub, local_field, out = sub )

    return field


def _shift_bool_array ( mask, dx, dy, dz, fill ):
    """ Returns an array the same shape as `mask` where result[i,j,k] =
    mask[i+dx, j+dy, k+dz] when that index is in bounds, else `fill`. Pure
    numpy slicing (no scipy.ndimage.shift), used by
    _bounded_edt_from_true_to_false() below. """
    out = np.full_like ( mask, fill )
    src_slices = [ ]
    dst_slices = [ ]
    for d, n in zip ( ( dx, dy, dz ), mask.shape ):
        if d >= 0:
            src_slices.append ( slice ( d, n ) )
            dst_slices.append ( slice ( 0, n - d ) )
        else:
            src_slices.append ( slice ( 0, n + d ) )
            dst_slices.append ( slice ( -d, n ) )
    out[ tuple ( dst_slices ) ] = mask[ tuple ( src_slices ) ]
    return out


def _offsets_within_radius ( spacing, max_dist ):
    """ Every integer voxel offset (dx,dy,dz), excluding (0,0,0), whose
    real-world distance (dx,dy,dz)*spacing is <= max_dist, sorted nearest
    first -- the brute-force "structuring element" used by
    _bounded_edt_from_true_to_false() below. Cheap: for the probe radii
    (~1.4 A) and grid spacings (~0.8 A) this module actually uses, this is
    on the order of a few dozen offsets, not thousands. """
    r_vox = int ( np.ceil ( max_dist / spacing ) )
    axis = np.arange ( -r_vox, r_vox + 1 )
    dx, dy, dz = np.meshgrid ( axis, axis, axis, indexing = "ij" )
    dx = dx.ravel ( ); dy = dy.ravel ( ); dz = dz.ravel ( )
    dist = spacing * np.sqrt ( ( dx.astype ( np.float64 ) ) ** 2 + ( dy.astype ( np.float64 ) ) ** 2 + ( dz.astype ( np.float64 ) ) ** 2 )
    keep = ( dist > 0.0 ) & ( dist <= max_dist )
    order = np.argsort ( dist[ keep ] )
    dx, dy, dz, dist = dx[ keep ][ order ], dy[ keep ][ order ], dz[ keep ][ order ], dist[ keep ][ order ]
    return list ( zip ( dx.tolist ( ), dy.tolist ( ), dz.tolist ( ), dist.tolist ( ) ) )


def _bounded_edt_from_true_to_false ( mask, spacing, max_dist ):
    """ [EN] scipy-free replacement for
    `scipy.ndimage.distance_transform_edt(mask) `, TRUNCATED at
    `max_dist` -- see this module's own top-of-file "NO SCIPY" note for
    why truncation is exact (not an approximation) for this module's own
    two use cases (a boolean "> probe_radius" decision, and a smooth
    field only ever sampled by marching cubes within ~1-2 voxels of its
    own zero-crossing).

    For every voxel where `mask` is True, returns the distance to the
    NEAREST voxel where `mask` is False (capped at `max_dist`, so any
    true distance beyond that just reads back as `max_dist` -- exactly
    matching scipy's own convention that this is what the caller wants
    whenever it only cares about "at least max_dist away"). Voxels
    outside the grid are treated as False (mask's own natural boundary --
    this module's grid is already padded well beyond the largest sphere
    it needs to resolve, so this never affects a real zero-crossing).
    Voxels where `mask` is already False get 0.0 (scipy's own
    convention). """
    dist = np.where ( mask, max_dist, 0.0 ).astype ( np.float64 )
    not_mask = ~mask
    for dx, dy, dz, d in _offsets_within_radius ( spacing, max_dist ):
        neighbor_not_mask = _shift_bool_array ( not_mask, dx, dy, dz, fill = True )
        np.minimum ( dist, np.where ( neighbor_not_mask, d, max_dist ), out = dist )
    dist[ ~mask ] = 0.0
    return dist


def _ses_field_from_sas_field ( field_sas, probe_radius, spacing ):
    """ EDTSurf's own two-pass Euclidean-distance-transform construction of
    the Solvent-Excluded (Connolly) surface from the already-computed SAS
    field -- see this module's own top-of-file docstring for the full
    derivation. Returns a continuous signed field (negative inside the
    SES solid, positive outside, real zero-crossing at the SES boundary)
    suitable for marching cubes at isovalue 0. """
    sas_occ = field_sas <= 0.0
    if not np.any ( sas_occ ) or np.all ( sas_occ ):
        # degenerate (no atoms reached this grid at all, or the whole grid
        # is SAS-occupied e.g. an absurdly coarse/small box) -- nothing
        # sensible to carve, fall back to the SAS field itself.
        return field_sas

    # distance from each SAS-occupied voxel to the nearest NON-SAS voxel
    # (i.e. the nearest point the probe's own bulk region can freely
    # occupy) -- exactly "how close is solvent bulk", per voxel, capped
    # just above probe_radius (that's the only threshold the boolean
    # decision below actually needs).
    dist_to_exterior = _bounded_edt_from_true_to_false ( sas_occ, spacing, probe_radius * 1.05 )
    ses_occ = sas_occ & ( dist_to_exterior > probe_radius )

    if not np.any ( ses_occ ):
        return field_sas
    if np.all ( ses_occ ):
        return -np.ones_like ( field_sas )

    # smooth signed field around the SES boundary for marching cubes'
    # own interpolation -- only needs to be accurate within a voxel or
    # two of the zero-crossing, so a small truncation radius is enough
    # (and much cheaper than the probe_radius-sized one above).
    near_field_cutoff = 2.0 * spacing
    outside = _bounded_edt_from_true_to_false ( ~ses_occ, spacing, near_field_cutoff )
    inside  = _bounded_edt_from_true_to_false (  ses_occ, spacing, near_field_cutoff )
    return outside - inside


class _MarchingCubesCorrupted ( Exception ):
    """ Raised internally when pDynamo3's own marching-cubes output looks
    corrupted (see this module's own top-of-file note for the confirmed
    root cause -- the compiled routine's vertex-buffer growth path loses
    data once the mesh needs more vertices than its own initial buffer
    allocation). Caught by build_surface_mesh(), which retries at a
    coarser grid spacing; never meant to escape this module. """
    pass


def _nearest_atom_colors ( points, atom_coords, atom_colors ):
    """ [EN] scipy-free (no cKDTree, see this module's own top-of-file
    "NO SCIPY" note) nearest-atom colour lookup -- plain chunked brute-
    force distance search. `points` is normally pDynamo3's own (already
    welded/shared, see top-of-file note) vertex array, so this typically
    runs on a few thousand rows against a few thousand atoms -- small
    enough that exact brute force (chunked only to bound peak memory for
    unusually large meshes) is both fast and, unlike the k=12-nearest-
    neighbour approximation this replaces, always exactly correct. """
    n_points = points.shape[0]
    out = np.empty ( ( n_points, 3 ), dtype = np.float32 )
    for start in range ( 0, n_points, _COLOR_CHUNK_SIZE ):
        chunk = points[ start : start + _COLOR_CHUNK_SIZE ]
        # (chunk_n, n_atoms) pairwise distances -- one chunk at a time so
        # this never allocates a full (n_points, n_atoms) matrix at once.
        diff = chunk[ :, None, : ] - atom_coords[ None, :, : ]
        dist2 = np.einsum ( "ijk,ijk->ij", diff, diff )
        nearest = np.argmin ( dist2, axis = 1 )
        out[ start : start + chunk.shape[0] ] = atom_colors[ nearest ]
    return out


_GHOST_EDGE_SIZE_FACTOR = 8.0   # matches _compute_valid_polygon_mask's own default exactly


def _edge_length_ghost_mask ( polygons_np, vertices_np, size_factor = _GHOST_EDGE_SIZE_FACTOR ):
    """ [EN] Same geometric ghost-triangle filter as surface_analysis_
    window.py's own `_compute_valid_polygon_mask()` (discard any triangle
    whose longest edge is a wild outlier -- >`size_factor` times -- vs the
    mesh's own median edge length; see that function's own docstring for
    the full derivation/history of this bug), reimplemented here to take
    already-converted numpy arrays directly instead of raw pDynamo3
    Array2D objects. Avoids converting `polygons`/`vertices` from pDynamo3
    a SECOND time (profiling found this redundant conversion -- one call
    inside the original `_compute_valid_polygon_mask`, a second one right
    after in this function -- was a real, avoidable cost per surface
    generated). Returns (mask, n_discarded), same as the function this
    mirrors. """
    if polygons_np.shape[0] == 0:
        return np.zeros ( 0, dtype = bool ), 0
    tri_verts = vertices_np[ polygons_np ]
    e01 = np.linalg.norm ( tri_verts[ :, 1, : ] - tri_verts[ :, 0, : ], axis = 1 )
    e12 = np.linalg.norm ( tri_verts[ :, 2, : ] - tri_verts[ :, 1, : ], axis = 1 )
    e20 = np.linalg.norm ( tri_verts[ :, 0, : ] - tri_verts[ :, 2, : ], axis = 1 )
    max_edge = np.maximum ( np.maximum ( e01, e12 ), e20 )
    typical = np.median ( max_edge )
    if typical == 0.0:
        typical = 1e-9
    mask = max_edge <= ( size_factor * typical )
    return mask, int ( polygons_np.shape[0] - mask.sum ( ) )


def _vertex_normals_from_mesh ( vertices_np, polygons_np ):
    """ [EN] PERFORMANCE FIX (profiling found pDynamo3's own
    `PolygonalSurface.MakePolygonNormals()`/`MakeVertexNormalsFromPolygonalNormals()`
    -- pure per-triangle/per-vertex PYTHON loops, see pScientific/Surfaces/
    PolygonalSurface.py -- were ~40% of this module's total per-frame
    time, ~0.4s on a real ~32k-vertex crambin.pdb surface). Replaces both
    calls with a single vectorised numpy computation of the EXACT SAME
    thing (confirmed by reading that pDynamo3 source directly): per-face
    normal = normalize(cross(v1-v0, v2-v0)); per-vertex normal = sum of
    its adjacent (already-unit-length) face normals, renormalized --
    an unweighted average of unit face normals, not area-weighted, same
    as pDynamo3's own convention. Degenerate (zero-area) faces/orphan
    vertices fall back to (0,0,1), matching this codebase's own existing
    convention for the same edge case (see vismol.utils.mesh_smoothing's
    taubin_smooth()). """
    tri_verts = vertices_np[ polygons_np ]
    face_normals = np.cross ( tri_verts[ :, 1, : ] - tri_verts[ :, 0, : ], tri_verts[ :, 2, : ] - tri_verts[ :, 0, : ] )
    face_lengths = np.linalg.norm ( face_normals, axis = 1 )
    safe_face = face_lengths > 1e-12
    face_normals[ safe_face ] /= face_lengths[ safe_face, None ]
    face_normals[ ~safe_face ] = 0.0

    vertex_normals = np.zeros_like ( vertices_np )
    np.add.at ( vertex_normals, polygons_np[ :, 0 ], face_normals )
    np.add.at ( vertex_normals, polygons_np[ :, 1 ], face_normals )
    np.add.at ( vertex_normals, polygons_np[ :, 2 ], face_normals )

    vertex_lengths = np.linalg.norm ( vertex_normals, axis = 1 )
    safe_vertex = vertex_lengths > 1e-12
    vertex_normals[ safe_vertex ] /= vertex_lengths[ safe_vertex, None ]
    vertex_normals[ ~safe_vertex ] = np.array ( [ 0.0, 0.0, 1.0 ] )
    return vertex_normals


def _isosurface_to_flat_mesh ( surface, atom_coords, atom_colors ):
    """ Converts a pDynamo3 PolygonalSurface into the flat (vertices,
    colors, indexes, normals) float32/uint32 tuple SurfaceRepresentation
    already expects (see representations.py) -- same triangle-extraction
    approach this app's own surface_parser()/surface_parser_mep() already
    use (reusing their own already-debugged helper directly, see this
    module's own top-of-file note for why), but coloured per-vertex by
    NEAREST ATOM (the user's own explicit request -- reuse EasyHybrid's
    existing per-atom display colors, not a new CPK table) instead of a
    single tiled color, and with NO Bohr->Angstrom division (this
    module's grid is already in Angstrom).

    Raises `_MarchingCubesCorrupted` if too large a fraction of the
    output is discarded by the (proven, geometric, reused from the QC
    surface tool) ghost-triangle filter -- see build_surface_mesh()'s
    retry loop, which is what actually handles it. """
    # Local import: an already-debugged helper from the QC surface tool,
    # reused as-is rather than duplicated -- imported lazily so this pure-
    # numpy/pDynamo3 module never pulls in gi/Gtk at module-load time.
    from gui.windows.analysis.surface_analysis_window import _pdynamo_array_to_numpy

    polygons = surface.polygons
    if polygons.rows == 0:
        empty3 = np.zeros ( ( 0, 3 ), dtype = np.float32 )
        return empty3.reshape ( -1 ), empty3.reshape ( -1 ), np.zeros ( 0, dtype = np.uint32 ), empty3.reshape ( -1 )

    polygons_np_all = _pdynamo_array_to_numpy ( polygons, np.int64 )
    vertices_np     = _pdynamo_array_to_numpy ( surface.vertices, np.float64 )

    valid_mask, n_discarded = _edge_length_ghost_mask ( polygons_np_all, vertices_np )
    ghost_fraction = n_discarded / float ( polygons.rows )
    if ghost_fraction > _MAX_GHOST_TRIANGLE_FRACTION:
        raise _MarchingCubesCorrupted (
            "{:.1%} of triangles discarded by the ghost-triangle filter (likely the "
            "compiled marching-cubes routine's vertex-buffer growth bug -- see "
            "molecular_surface.py's own top-of-file note)".format ( ghost_fraction )
        )
    polygons_np = polygons_np_all[ valid_mask ]

    normals_np    = _vertex_normals_from_mesh ( vertices_np, polygons_np )
    vertex_colors = _nearest_atom_colors ( vertices_np, atom_coords, atom_colors )

    tri_verts  = vertices_np[ polygons_np ]      # (n_valid_tri, 3, 3)
    tri_norms  = normals_np[ polygons_np ]
    tri_colors = vertex_colors[ polygons_np ]

    vertices_out = tri_verts.reshape ( -1 ).astype ( np.float32 )
    normals_out  = tri_norms.reshape ( -1 ).astype ( np.float32 )
    colors_out   = tri_colors.reshape ( -1 ).astype ( np.float32 )

    indexes_out = np.arange ( vertices_out.shape[0] // 3, dtype = np.uint32 )
    return vertices_out, colors_out, indexes_out, normals_out


def _build_surface_mesh_from_atoms ( coords, vdw_radii, colors, surface_type,
                                      probe_radius = DEFAULT_PROBE_RADIUS,
                                      grid_spacing = DEFAULT_GRID_SPACING,
                                      taubin_iterations = 0 ):
    """ [EN] The actual marching-cubes pipeline, operating ONLY on plain
    numpy arrays (one frame's worth of atom coordinates/radii/colors) --
    no VismolObject/GTK dependency at all. Split out from build_surface_
    mesh() (below, now a thin per-object/per-frame wrapper around this)
    so it can ALSO be called as a multiprocessing worker (see
    _generate_surface_frame_job()/build_surface_trajectory() below) --
    plain numpy arrays are trivially picklable across the process
    boundary, unlike a live VismolObject/pDynamo3 C object.

    `surface_type`: "vws" (Van der Waals), "sas" (Solvent-Accessible), or
    "ses" (Solvent-Excluded / Connolly) -- see this module's own top-of-
    file docstring for what each one means and how it's computed.
    `probe_radius`: only used for "sas"/"ses" (ignored for "vws").
    `taubin_iterations`: 0 (default) leaves the geometry exactly as
    marching cubes produced it (still smoothly NORMAL-shaded, see this
    module's own top-of-file note -- that part needs no smoothing pass at
    all). >0 additionally runs `vismol.utils.mesh_smoothing`'s own
    already-proven Taubin (lambda|mu) implementation for that many
    iterations -- the user's own explicit request, "usar o algoritmo de
    Taubin (15 iteracoes) para suavizar a superficie" -- reused as-is
    rather than reimplemented (it already does exactly this: weld
    coincident vertex slots, smooth positions without the shrinkage plain
    Laplacian smoothing causes, then recompute per-vertex normals from
    the smoothed geometry).

    [EN] `grid_spacing` is a STARTING point, not a guarantee -- see this
    module's own top-of-file note on the real, confirmed pDynamo3
    marching-cubes vertex-buffer corruption bug. If the requested spacing
    would trigger it, this automatically retries at a coarser spacing
    (up to `_MAX_RESOLUTION_RETRIES` times) until the output is clean. """
    surface_type = surface_type.lower ( )
    if surface_type not in ( "vws", "sas", "ses" ):
        raise ValueError ( "unknown surface_type {!r} -- expected 'vws', 'sas' or 'ses'".format ( surface_type ) )

    if len ( coords ) == 0:
        empty3 = np.zeros ( ( 0, 3 ), dtype = np.float32 )
        return empty3.reshape ( -1 ), empty3.reshape ( -1 ), np.zeros ( 0, dtype = np.uint32 ), empty3.reshape ( -1 )

    spacing = grid_spacing
    last_error = None
    mesh = None
    for attempt in range ( _MAX_RESOLUTION_RETRIES + 1 ):
        if surface_type == "vws":
            max_radius = vdw_radii.max ( )
            lower, dims = _grid_bounds ( coords, max_radius, spacing )
            field = _signed_union_of_spheres_field ( coords, vdw_radii, lower, dims, spacing )
        else:
            sas_radii = vdw_radii + probe_radius
            max_radius = sas_radii.max ( )
            lower, dims = _grid_bounds ( coords, max_radius, spacing )
            field_sas = _signed_union_of_spheres_field ( coords, sas_radii, lower, dims, spacing )
            if surface_type == "sas":
                field = field_sas
            else:   # "ses"
                field = _ses_field_from_sas_field ( field_sas, probe_radius, spacing )

        nx, ny, nz = dims
        grid = RegularGrid.FromDimensionData ( [
            { "bins": nx, "binSize": spacing, "lower": float ( lower[0] ) },
            { "bins": ny, "binSize": spacing, "lower": float ( lower[1] ) },
            { "bins": nz, "binSize": spacing, "lower": float ( lower[2] ) },
        ] )
        flat_array = Array.FromIterable ( field.reshape ( -1 ).tolist ( ) )
        dataND     = Reshape ( flat_array, ( nx, ny, nz ), resultClass = RealArrayND )

        surface = _isosurface_from_grid ( grid, dataND, 0.0 )

        try:
            mesh = _isosurface_to_flat_mesh ( surface, coords, colors )
            break
        except _MarchingCubesCorrupted as error:
            last_error = error
            spacing *= _GRID_SPACING_RETRY_GROWTH
    else:
        raise RuntimeError (
            "molecular surface generation failed after {} attempts (last: {}) -- "
            "pDynamo3's marching-cubes routine kept exceeding its own vertex-buffer "
            "corruption threshold even at a coarse grid spacing of {:.2f} A; try a "
            "smaller atom selection.".format ( _MAX_RESOLUTION_RETRIES + 1, last_error, spacing )
        )

    if taubin_iterations and taubin_iterations > 0:
        from vismol.utils import mesh_smoothing
        vertices_out, colors_out, indexes_out, normals_out = mesh
        smoothed_vertices, smoothed_normals = mesh_smoothing.taubin_smooth (
            np.ascontiguousarray ( vertices_out, dtype = np.float32 ),
            np.ascontiguousarray ( indexes_out,  dtype = np.uint32  ),
            iterations = int ( taubin_iterations ),
        )
        mesh = ( smoothed_vertices, colors_out, indexes_out, smoothed_normals )

    return mesh


def _triangle_areas ( vertices ):
    """ [EN] (n_tri,) array of triangle areas (same units as the input
    coordinates -- Angstrom^2 for every caller in this module), from a
    flat (3*3*n_tri,) vertex array -- this module's own (and
    SurfaceRepresentation's own) "triangle soup" layout: 3 OWN vertices
    per triangle, no shared indexing (see this module's own top-of-file
    note on why the final OUTPUT is always in this format, even though
    pDynamo3's own vertices/polygons are a real welded mesh internally).
    Area = 0.5 * |cross(v1-v0, v2-v0)|, exact for triangles. """
    tri = vertices.reshape ( -1, 3, 3 )
    cross = np.cross ( tri[ :, 1, : ] - tri[ :, 0, : ], tri[ :, 2, : ] - tri[ :, 0, : ] )
    return 0.5 * np.linalg.norm ( cross, axis = 1 )


def compute_surface_area ( mesh, triangle_mask = None ):
    """ [EN] Total area (Angstrom^2) of a (vertices, colors, indexes,
    normals) mesh tuple -- the exact structure build_surface_mesh()/
    build_surface_trajectory() already return, and what `vismol_object.
    surface_trajectory[frame][surf_name]` already stores -- so this can
    be called directly on any already-generated surface, no regeneration
    needed. Simple sum of per-triangle areas (see _triangle_areas()
    above); exact (not an approximation) for this module's own always-
    triangulated marching-cubes output.

    `triangle_mask`: optional (n_tri,) boolean array -- when given, only
    the SELECTED triangles are summed, instead of the whole surface. Not
    used by any caller yet -- the user's own explicit request, 2026-09-24:
    "ja prepare para o calculo fazer para um pedaco dela tambem" (a
    future partial/patch area query, e.g. restricted to a specific
    residue/selection/nearest-atom subset) -- but the parameter already
    exists so that can be added later as a pure caller-side change (build
    the right boolean mask, e.g. from a per-triangle nearest-atom lookup
    -- see _nearest_atom_colors()'s own nearest-atom index, computed but
    currently discarded after colouring -- and pass it here), with ZERO
    changes needed to this function's own signature or the "whole
    surface" call sites that already exist. """
    areas = _triangle_areas ( mesh[0] )
    if triangle_mask is not None:
        areas = areas[ triangle_mask ]
    return float ( areas.sum ( ) )


def build_surface_mesh ( vismol_object, surface_type = "ses",
                          probe_radius = DEFAULT_PROBE_RADIUS,
                          grid_spacing = DEFAULT_GRID_SPACING,
                          frame = None, taubin_iterations = 0, color_mode = "atom" ):
    """ Single-frame convenience wrapper -- gathers atom coordinates/radii/
    colors for ONE frame of `vismol_object` (`frame=None` resolves to
    "whichever frame is currently active", same convention every other
    representation in this app already uses) and builds its
    (vertices, colors, indexes, normals) flat mesh tuple via
    _build_surface_mesh_from_atoms() above (see its own docstring for what
    every parameter means). `color_mode`: see _resolve_atom_colors()'s own
    docstring ("atom"/"chain"/"carbon"). For a WHOLE trajectory (every
    frame, not just one), use build_surface_trajectory() below instead. """
    coords, vdw_radii, colors = _gather_atoms ( vismol_object, frame = frame, color_mode = color_mode )
    return _build_surface_mesh_from_atoms (
        coords, vdw_radii, colors, surface_type,
        probe_radius = probe_radius, grid_spacing = grid_spacing, taubin_iterations = taubin_iterations,
    )


def _generate_surface_frame_job ( job ):
    """ [EN] Top-level (module-scope, hence picklable) multiprocessing
    worker -- builds ONE trajectory frame's mesh. `job` is a plain tuple
    (not a dict/object) of picklable values: (coords, vdw_radii, colors,
    surface_type, probe_radius, grid_spacing, taubin_iterations). Used by
    build_surface_trajectory() below via `multiprocessing.Pool.map()` --
    the exact same parallelisation pattern surface_analysis_window.py's
    own generate_grid_parallel() already uses for QC surface trajectories
    (see that file, e.g. its density branch: `p = multiprocessing.Pool
    (processes=multiprocessing.cpu_count()); results = p.map
    (generate_grid_parallel, joblist)`). """
    coords, vdw_radii, colors, surface_type, probe_radius, grid_spacing, taubin_iterations = job
    return _build_surface_mesh_from_atoms (
        coords, vdw_radii, colors, surface_type,
        probe_radius = probe_radius, grid_spacing = grid_spacing, taubin_iterations = taubin_iterations,
    )


def build_surface_trajectory ( vismol_object, surface_type = "ses",
                                probe_radius = DEFAULT_PROBE_RADIUS,
                                grid_spacing = DEFAULT_GRID_SPACING,
                                taubin_iterations = 0, surf_name = "surface", color_mode = "atom" ):
    """ [EN] Builds a REAL mesh for EVERY frame of `vismol_object`'s own
    trajectory (`vismol_object.frames`), not a single static frame copied
    into every slot -- the user's own explicit follow-up request,
    2026-09-23: "extenda para todos os frames". Atom radii/colors are
    gathered once (they don't vary by frame); only coordinates are
    re-gathered per frame. Runs sequentially for a single-frame object
    (the common case -- most loaded structures have no trajectory, so
    this stays exactly as fast as before), or in PARALLEL via
    `multiprocessing.Pool` when there is more than one frame -- see
    _generate_surface_frame_job()'s own docstring for why this mirrors
    the QC surface tool's own established pattern, not a new one.
    `color_mode`: see _resolve_atom_colors()'s own docstring
    ("atom"/"chain"/"carbon").

    Returns a list of `{surf_name: (vertices, colors, indexes, normals)}`
    dicts, one per frame, in frame order -- ready to assign directly to a
    VismolObject's own `.surface_trajectory` (the exact structure
    SurfaceRepresentation.draw_representation() already reads per frame,
    see representations.py). """
    atom_ids = sorted ( vismol_object.atoms.keys ( ) )
    atoms = [ vismol_object.atoms[i] for i in atom_ids ]
    n_frames = int ( vismol_object.frames.shape[0] )
    if not atoms:
        return [ { } for _ in range ( n_frames ) ]

    periodic_table = vismol_object.vm_session.periodic_table
    vdw_radii = np.array (
        [ periodic_table.get_vdw_radius ( symbol = atom.symbol ) for atom in atoms ],
        dtype = np.float64,
    )
    colors = _resolve_atom_colors ( vismol_object, atom_ids, atoms, color_mode )

    joblist = [ ]
    for frame in range ( n_frames ):
        coords = np.array ( [ atom.coords ( frame ) for atom in atoms ], dtype = np.float64 )
        joblist.append ( ( coords, vdw_radii, colors, surface_type, probe_radius, grid_spacing, taubin_iterations ) )

    if n_frames <= 1:
        meshes = [ _generate_surface_frame_job ( joblist[0] ) ] if joblist else [ ]
    else:
        import multiprocessing
        with multiprocessing.Pool ( processes = min ( multiprocessing.cpu_count ( ), n_frames ) ) as pool:
            meshes = pool.map ( _generate_surface_frame_job, joblist )

    return [ { surf_name: mesh } for mesh in meshes ]
