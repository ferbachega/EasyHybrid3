#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  cavities.py
#
#  Cavity and pocket detection on a grid (2026-10-30), the method of
#  KVFinder (Oliveira et al., BMC Bioinformatics 2014; parKVFinder,
#  Guerra et al., SoftwareX 2020), in numpy only (no scipy, as
#  molecular_surface.py) and with no GTK, so it runs in worker processes:
#
#    - Probe In (1.4 A): the points a water-sized probe can reach, the
#      union of probe spheres whose centres touch no atom (vdW + probe).
#      Buried voids where the probe fits are included.
#    - Probe Out (4 A): the same with a large probe rolled from OUTSIDE
#      the molecule (only centres connected to the border of the grid):
#      the solvent bulk. The molecule's "hull" is what it cannot reach.
#    - cavity points = reachable by Probe In, not by Probe Out, and more
#      than the removal distance (2.4 A) from the bulk (shallow grooves of
#      the surface are not cavities);
#    - cavities = connected groups of points (6-neighbours); groups smaller
#      than the minimum volume are dropped. Volume = points x spacing^3.
#
#  Per cavity: volume, area (marching-cubes mesh of the cavity), depth
#  (distance of its points to the mouth: the points next to the bulk; a
#  cavity with no mouth is buried) and the atoms lining it.
#
#  Along a trajectory a pocket is followed by its lining atoms (atom
#  indices do not change between frames; coordinates may drift):
#  pocket_volume_trajectory().

import numpy as np

from util.molecular_surface import ( _signed_union_of_spheres_field, _shift_bool_array,
                                     _bounded_edt_from_true_to_false, _isosurface_from_grid,
                                     _isosurface_to_flat_mesh, compute_surface_area,
                                     _MarchingCubesCorrupted )

DEFAULT_SPACING          = 0.6    # A, KVFinder's default step
DEFAULT_PROBE_IN         = 1.4    # A
DEFAULT_PROBE_OUT        = 4.0    # A
DEFAULT_REMOVAL_DISTANCE = 2.4    # A
DEFAULT_MIN_VOLUME       = 5.0    # A^3
LINING_MARGIN            = 0.6    # A beyond the vdW radius: atoms touching the cavity points

# . Bondi (J. Phys. Chem. 1964) van der Waals radii, A, the usual ones for
#   molecular volumes; the display radii of EasyHybrid are larger (C 1.95)
#   and would shrink every cavity. Others: the caller's fallback.
BONDI_RADII = { "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80, "Cl": 1.75,
                "Br": 1.85, "I": 1.98, "Se": 1.90, "Si": 2.10, "Na": 2.27, "K": 2.75, "Mg": 1.73,
                "Zn": 1.39, "Cu": 1.40, "Ni": 1.63, "Fe": 1.94, "Mn": 1.97, "Co": 1.92, "Ca": 2.31 }


def radii_for ( symbols, fallback = None ):
    """ Bondi radii of the element symbols; `fallback(symbol)` (or 1.8 A) for others. """
    out = [ ]
    for s in symbols:
        s = s.capitalize ( ) if s else s
        if s in BONDI_RADII:
            out.append ( BONDI_RADII[ s ] )
        else:
            out.append ( float ( fallback ( s ) ) if fallback else 1.8 )
    return np.array ( out, dtype = np.float64 )


_SIX = ( ( 1, 0, 0 ), ( -1, 0, 0 ), ( 0, 1, 0 ), ( 0, -1, 0 ), ( 0, 0, 1 ), ( 0, 0, -1 ) )

CAVITY_COLORS = [ ( 0.95, 0.35, 0.25 ), ( 0.25, 0.55, 0.95 ), ( 0.30, 0.80, 0.35 ), ( 0.95, 0.75, 0.20 ),
                  ( 0.70, 0.40, 0.90 ), ( 0.25, 0.80, 0.80 ), ( 0.95, 0.45, 0.70 ), ( 0.60, 0.45, 0.30 ) ]


# --------------------------------------------------------------------------- grid morphology
def _ball_offsets ( radius_vox ):
    """ Integer offsets (dx,dy,dz) with |offset| <= radius_vox, (0,0,0) excluded. """
    r = int ( np.floor ( radius_vox ) )
    axis = np.arange ( -r, r + 1 )
    dx, dy, dz = np.meshgrid ( axis, axis, axis, indexing = "ij" )
    keep = ( dx ** 2 + dy ** 2 + dz ** 2 <= radius_vox ** 2 + 1e-9 ) & ( ( dx != 0 ) | ( dy != 0 ) | ( dz != 0 ) )
    return list ( zip ( dx[ keep ].tolist ( ), dy[ keep ].tolist ( ), dz[ keep ].tolist ( ) ) )


def dilate_ball ( mask, radius, spacing ):
    """ Points within `radius` (A) of a True point. Exact ball for small
        radii; for large ones, successive dilations by smaller balls (the
        Minkowski sum of balls is a ball), within about one grid step. """
    if radius <= 0.0:
        return mask.copy ( )
    radius_vox = radius / spacing
    steps = max ( 1, int ( np.ceil ( radius_vox / 4.0 ) ) )
    out = mask.copy ( )
    for _ in range ( steps ):
        src = out.copy ( )
        for dx, dy, dz in _ball_offsets ( radius_vox / steps ):
            out |= _shift_bool_array ( src, dx, dy, dz, fill = False )
    return out


def _border_connected ( mask ):
    """ The part of `mask` connected (6-neighbours) to the border of the grid. """
    seed = np.zeros_like ( mask )
    seed[ 0, :, : ] = seed[ -1, :, : ] = True
    seed[ :, 0, : ] = seed[ :, -1, : ] = True
    seed[ :, :, 0 ] = seed[ :, :, -1 ] = True
    reached = seed & mask
    while True:
        grown = reached.copy ( )
        for d in _SIX:
            grown |= _shift_bool_array ( reached, *d, fill = False )
        grown &= mask
        if grown.sum ( ) == reached.sum ( ):
            return reached
        reached = grown


def label_components ( mask ):
    """ (labels, n): connected components of `mask` (6-neighbours), labels
        1..n in the grid (0 elsewhere). Union-find on the True points only. """
    idx = np.flatnonzero ( mask )
    if idx.size == 0:
        return np.zeros ( mask.shape, dtype = np.int32 ), 0
    nx, ny, nz = mask.shape
    position = -np.ones ( mask.size, dtype = np.int64 )
    position[ idx ] = np.arange ( idx.size )
    i, j, k = np.unravel_index ( idx, mask.shape )
    edges_a, edges_b = [ ], [ ]
    for axis, ( a, n ) in enumerate ( ( ( i, nx ), ( j, ny ), ( k, nz ) ) ):
        ok = a + 1 < n
        step = ( ny * nz, nz, 1 )[ axis ]
        neighbour = idx[ ok ] + step
        linked = mask.reshape ( -1 )[ neighbour ]
        edges_a.append ( position[ idx[ ok ][ linked ] ] )
        edges_b.append ( position[ neighbour[ linked ] ] )
    ea = np.concatenate ( edges_a )
    eb = np.concatenate ( edges_b )
    parent = np.arange ( idx.size )
    while True:
        pa, pb = parent[ ea ], parent[ eb ]
        low = np.minimum ( pa, pb )
        new = parent.copy ( )
        np.minimum.at ( new, pa, low )
        np.minimum.at ( new, pb, low )
        new = new[ new ]                       # pointer jumping
        while True:
            jumped = new[ new ]
            if np.array_equal ( jumped, new ):
                break
            new = jumped
        if np.array_equal ( new, parent ):
            break
        parent = new
    roots, compact = np.unique ( parent, return_inverse = True )
    labels = np.zeros ( mask.size, dtype = np.int32 )
    labels[ idx ] = compact + 1
    return labels.reshape ( mask.shape ), int ( roots.size )


# --------------------------------------------------------------------------- detection
def _grid ( coords, radii, spacing, probe_out ):
    pad = float ( radii.max ( ) ) + 2.0 * probe_out + 3.0 * spacing
    lower = coords.min ( axis = 0 ) - pad
    upper = coords.max ( axis = 0 ) + pad
    dims = tuple ( int ( d ) for d in np.ceil ( ( upper - lower ) / spacing ).astype ( np.int64 ) + 1 )
    return lower, dims


def grid_points ( result, flat_indices ):
    """ Coordinates (A) of grid points given by flat indices. """
    i, j, k = np.unravel_index ( flat_indices, result[ "dims" ] )
    return result[ "lower" ] + 0.5 * result[ "spacing" ] + np.stack ( ( i, j, k ), axis = 1 ) * result[ "spacing" ]


def detect_cavities ( coords, radii, spacing = DEFAULT_SPACING, probe_in = DEFAULT_PROBE_IN,
                      probe_out = DEFAULT_PROBE_OUT, removal_distance = DEFAULT_REMOVAL_DISTANCE,
                      min_volume = DEFAULT_MIN_VOLUME, box = None, lining = True ):
    """ Cavities of the atoms (coords (n,3) A, vdW radii (n,) A).

        box: optional (centre (3,), half sizes (3,)) in A: only points inside
        are kept (a pocket of interest); the bulk is still computed on the
        whole molecule.

        Returns {"lower", "dims", "spacing", "labels" (grid int32, 0 none,
        1..n by decreasing volume), "cavities": [ {id, points, volume, area
        (None until a mesh is made), depth_max, depth_mean, buried, centre,
        lining_atoms} ], parameters}. """
    coords = np.asarray ( coords, dtype = np.float64 )
    radii = np.asarray ( radii, dtype = np.float64 )
    lower, dims = _grid ( coords, radii, spacing, probe_out )

    # . Probe In: centres touching no atom, then the spheres around them
    centres_in = _signed_union_of_spheres_field ( coords, radii + probe_in, lower, dims, spacing ) > 0.0
    reach_in = dilate_ball ( centres_in, probe_in, spacing )
    # . Probe Out from outside only: the bulk solvent
    centres_out = _signed_union_of_spheres_field ( coords, radii + probe_out, lower, dims, spacing ) > 0.0
    bulk = dilate_ball ( _border_connected ( centres_out ), probe_out, spacing )

    cavity = reach_in & ~bulk
    near_bulk = dilate_ball ( bulk, removal_distance, spacing )
    cavity &= ~near_bulk
    if box is not None:
        centre, half = ( np.asarray ( v, dtype = np.float64 ) for v in box )
        inside = [ np.abs ( lower[a] + 0.5 * spacing + np.arange ( dims[a] ) * spacing - centre[a] ) <= half[a] for a in range ( 3 ) ]
        cavity &= inside[0][ :, None, None ] & inside[1][ None, :, None ] & inside[2][ None, None, : ]

    labels, n = label_components ( cavity )
    # . the mouth: cavity points next to the removed shell or the bulk (not next to atoms)
    open_space = reach_in & near_bulk
    next_to_open = np.zeros_like ( cavity )
    for d in _SIX:
        next_to_open |= _shift_bool_array ( open_space, *d, fill = False )
    mouth = cavity & next_to_open

    flat_labels = labels.reshape ( -1 )
    order = np.argsort ( flat_labels[ flat_labels > 0 ], kind = "stable" )
    indices = np.flatnonzero ( flat_labels > 0 )[ order ]
    counts = np.bincount ( flat_labels, minlength = n + 1 )[ 1: ]
    groups = np.split ( indices, np.cumsum ( counts )[ :-1 ] ) if n else [ ]
    mouth_flat = mouth.reshape ( -1 )

    voxel = spacing ** 3
    cavities = [ ]
    for points in groups:
        volume = points.size * voxel
        if volume < min_volume:
            continue
        xyz = grid_points ( { "lower": lower, "dims": dims, "spacing": spacing }, points )
        mouth_xyz = xyz[ mouth_flat[ points ] ]
        if len ( mouth_xyz ):
            depth = _min_distances ( xyz, mouth_xyz )
            depth_max, depth_mean, buried = float ( depth.max ( ) ), float ( depth.mean ( ) ), False
        else:
            depth_max = depth_mean = 0.0
            buried = True
        cavities.append ( { "points": points, "volume": float ( volume ), "area": None,
                            "depth_max": depth_max, "depth_mean": depth_mean, "buried": buried,
                            "centre": xyz.mean ( axis = 0 ), "lining_atoms": None } )

    cavities.sort ( key = lambda c: -c[ "volume" ] )
    new_labels = np.zeros ( labels.size, dtype = np.int32 )
    for number, c in enumerate ( cavities, start = 1 ):
        c[ "id" ] = number
        new_labels[ c[ "points" ] ] = number
    result = { "lower": lower, "dims": dims, "spacing": spacing, "labels": new_labels.reshape ( dims ),
               "cavities": cavities,
               "parameters": { "spacing": spacing, "probe_in": probe_in, "probe_out": probe_out,
                               "removal_distance": removal_distance, "min_volume": min_volume } }
    if lining:
        assign_lining_atoms ( result, coords, radii )
    return result


def _min_distances ( points, targets, chunk = 2048 ):
    """ For each point, the distance to the nearest target (brute force, in chunks). """
    out = np.empty ( len ( points ) )
    for start in range ( 0, len ( points ), chunk ):
        block = points[ start:start + chunk ]
        best = np.full ( len ( block ), np.inf )
        for t0 in range ( 0, len ( targets ), chunk ):
            d2 = ( ( block[ :, None, : ] - targets[ None, t0:t0 + chunk, : ] ) ** 2 ).sum ( axis = 2 )
            np.minimum ( best, d2.min ( axis = 1 ), out = best )
        out[ start:start + chunk ] = np.sqrt ( best )
    return out


def assign_lining_atoms ( result, coords, radii, margin = LINING_MARGIN ):
    """ cavity["lining_atoms"]: the atoms whose vdW sphere, grown by
        `margin`, contains a point of the cavity (sorted atom indices). """
    labels, lower, spacing, dims = result[ "labels" ], result[ "lower" ], result[ "spacing" ], result[ "dims" ]
    lining = { c[ "id" ]: [ ] for c in result[ "cavities" ] }
    dims_arr = np.array ( dims )
    for atom, ( centre, radius ) in enumerate ( zip ( coords, radii ) ):
        reach = radius + margin
        lo = np.maximum ( np.floor ( ( centre - reach - lower ) / spacing - 0.5 ).astype ( np.int64 ), 0 )
        hi = np.minimum ( np.ceil ( ( centre + reach - lower ) / spacing - 0.5 ).astype ( np.int64 ) + 1, dims_arr )
        if np.any ( hi <= lo ):
            continue
        block = labels[ lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2] ]
        if not block.any ( ):
            continue
        axes = [ lower[a] + 0.5 * spacing + np.arange ( lo[a], hi[a] ) * spacing - centre[a] for a in range ( 3 ) ]
        d2 = axes[0][ :, None, None ] ** 2 + axes[1][ None, :, None ] ** 2 + axes[2][ None, None, : ] ** 2
        for number in np.unique ( block[ ( block > 0 ) & ( d2 <= reach * reach ) ] ):
            lining[ int ( number ) ].append ( atom )
    for c in result[ "cavities" ]:
        c[ "lining_atoms" ] = np.array ( lining[ c[ "id" ] ], dtype = np.int64 )
    return result


# --------------------------------------------------------------------------- meshes
def cavity_mesh ( result, cavity_id, color = None, padding = 3 ):
    """ (vertices, colors, indexes, normals) mesh of one cavity (the format
        of SurfaceRepresentation), its area set in the cavity record. The
        surface passes half a grid step outside the cavity points. """
    from pScientific.Arrays import Array, Reshape, RealArrayND
    from pScientific.Geometry3 import RegularGrid

    cavity = result[ "cavities" ][ cavity_id - 1 ]
    spacing, lower = result[ "spacing" ], result[ "lower" ]
    ijk = np.stack ( np.unravel_index ( cavity[ "points" ], result[ "dims" ] ), axis = 1 )
    lo = ijk.min ( axis = 0 ) - padding
    hi = ijk.max ( axis = 0 ) + padding + 1
    dims = tuple ( int ( v ) for v in hi - lo )
    inside = np.zeros ( dims, dtype = bool )
    inside[ tuple ( ( ijk - lo ).T ) ] = True
    cutoff = 2.0 * spacing
    # . -spacing on the last inside point, +spacing on the first outside one: the zero
    #   level is half a step outside the points, the faces of their cells (volume = n h^3)
    field = _bounded_edt_from_true_to_false ( ~inside, spacing, cutoff ) - _bounded_edt_from_true_to_false ( inside, spacing, cutoff )
    sub_lower = lower + lo * spacing
    grid = RegularGrid.FromDimensionData ( [ { "bins": dims[a], "binSize": spacing, "lower": float ( sub_lower[a] ) }
                                             for a in range ( 3 ) ] )
    data = Reshape ( Array.FromIterable ( field.reshape ( -1 ).tolist ( ) ), dims, resultClass = RealArrayND )
    rgb = np.array ( [ color or CAVITY_COLORS[ ( cavity_id - 1 ) % len ( CAVITY_COLORS ) ] ], dtype = np.float32 )
    surface = _isosurface_from_grid ( grid, data, 0.0 )
    try:
        mesh = _isosurface_to_flat_mesh ( surface, cavity[ "centre" ][ None, : ], rgb )
    except _MarchingCubesCorrupted:
        mesh = None
    if mesh is not None:
        cavity[ "area" ] = compute_surface_area ( mesh )
    return mesh


def points_mesh ( result, cavity_id, radius = 0.22, color = None ):
    """ The cavity points as small spheres (octahedra subdivided once), in the
        mesh format of cavity_mesh(). """
    xyz = grid_points ( result, result[ "cavities" ][ cavity_id - 1 ][ "points" ] ).astype ( np.float32 )
    v = np.array ( [ [ 1, 0, 0 ], [ -1, 0, 0 ], [ 0, 1, 0 ], [ 0, -1, 0 ], [ 0, 0, 1 ], [ 0, 0, -1 ] ], dtype = np.float32 )
    faces = [ ( 0, 2, 4 ), ( 2, 1, 4 ), ( 1, 3, 4 ), ( 3, 0, 4 ), ( 2, 0, 5 ), ( 1, 2, 5 ), ( 3, 1, 5 ), ( 0, 3, 5 ) ]
    tris = [ ]
    for a, b, c in faces:                      # . one subdivision: 32 triangles
        A, B, C = v[ a ], v[ b ], v[ c ]
        ab, bc, ca = ( A + B ), ( B + C ), ( C + A )
        ab, bc, ca = ( x / np.linalg.norm ( x ) for x in ( ab, bc, ca ) )
        tris += [ ( A, ab, ca ), ( ab, B, bc ), ( ca, bc, C ), ( ab, bc, ca ) ]
    unit = np.array ( tris, dtype = np.float32 )                       # (32, 3, 3) on the unit sphere
    verts = ( xyz[ :, None, None, : ] + radius * unit[ None ] ).reshape ( -1, 3 )
    normals = np.broadcast_to ( unit[ None ], ( len ( xyz ), ) + unit.shape ).reshape ( -1, 3 )
    rgb = np.asarray ( color or CAVITY_COLORS[ ( cavity_id - 1 ) % len ( CAVITY_COLORS ) ], dtype = np.float32 )
    colors = np.broadcast_to ( rgb, verts.shape )
    return ( verts.reshape ( -1 ).astype ( np.float32 ), np.ascontiguousarray ( colors ).reshape ( -1 ).astype ( np.float32 ),
             np.arange ( len ( verts ), dtype = np.uint32 ), np.ascontiguousarray ( normals ).reshape ( -1 ).astype ( np.float32 ) )


def box_around ( coords, margin ):
    """ (centre, half sizes) of the box around the coordinates, grown by margin (A). """
    coords = np.asarray ( coords, dtype = np.float64 )
    lo, hi = coords.min ( axis = 0 ) - margin, coords.max ( axis = 0 ) + margin
    return ( lo + hi ) / 2.0, ( hi - lo ) / 2.0


# --------------------------------------------------------------------------- reports and files
def lining_residues ( cavity, atom_labels ):
    """ Distinct residues of the lining atoms, in order. atom_labels: list of
        (chain, residue name, residue number) per atom index. """
    seen, out = set ( ), [ ]
    for atom in cavity[ "lining_atoms" ]:
        key = atom_labels[ atom ]
        if key not in seen:
            seen.add ( key )
            out.append ( key )
    return sorted ( out, key = lambda r: ( str ( r[0] ), r[2] ) )


def write_points_pdb ( result, path, cavity_ids = None ):
    """ Cavity points as HETATM records, residue name CAV, residue number =
        cavity id (one per cavity; all cavities when cavity_ids is None). """
    lines, serial = [ ], 1
    for c in result[ "cavities" ]:
        if cavity_ids is not None and c[ "id" ] not in cavity_ids:
            continue
        for x, y, z in grid_points ( result, c[ "points" ] ):
            lines.append ( "HETATM{:5d}  K   CAV Z{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}           K".format (
                serial % 100000, c[ "id" ] % 10000, x, y, z, 0.0 ) )
            serial += 1
    lines.append ( "END" )
    with open ( path, "w" ) as f:
        f.write ( "\n".join ( lines ) + "\n" )
    return path


def write_table_csv ( result, path, atom_labels = None ):
    import csv
    with open ( path, "w", newline = "" ) as f:
        w = csv.writer ( f )
        w.writerow ( [ "cavity", "volume_A3", "area_A2", "depth_max_A", "depth_mean_A", "buried", "centre_x",
                       "centre_y", "centre_z", "lining_residues" ] )
        for c in result[ "cavities" ]:
            residues = ""
            if atom_labels is not None:
                residues = " ".join ( "{}:{}{}".format ( ch, name, num ) for ch, name, num in lining_residues ( c, atom_labels ) )
            w.writerow ( [ c[ "id" ], "{:.1f}".format ( c[ "volume" ] ),
                           "" if c[ "area" ] is None else "{:.1f}".format ( c[ "area" ] ),
                           "{:.2f}".format ( c[ "depth_max" ] ), "{:.2f}".format ( c[ "depth_mean" ] ),
                           "yes" if c[ "buried" ] else "no" ] + [ "{:.3f}".format ( v ) for v in c[ "centre" ] ] + [ residues ] )
    return path


# --------------------------------------------------------------------------- trajectories
def match_pocket ( result, reference_atoms, min_overlap = 0.3 ):
    """ The cavities of `result` that are the reference pocket: those whose
        lining atoms overlap the reference lining atoms (fraction of the
        cavity's own lining in the reference set >= min_overlap). Returns
        (volume summed over them, [ids], overlap of the best one). """
    reference = set ( int ( a ) for a in reference_atoms )
    ids, volume, best = [ ], 0.0, 0.0
    for c in result[ "cavities" ]:
        lining = set ( c[ "lining_atoms" ].tolist ( ) )
        if not lining:
            continue
        shared = len ( lining & reference )
        overlap = shared / float ( min ( len ( lining ), len ( reference ) ) )
        best = max ( best, overlap )
        if overlap >= min_overlap:
            ids.append ( c[ "id" ] )
            volume += c[ "volume" ]
    return volume, ids, best


def _pocket_frame_job ( job ):
    coords, radii, options, reference_atoms, min_overlap = job
    options = dict ( options )
    box_atoms = options.pop ( "box_atoms", None )        # . (indices, margin): a box that follows those atoms
    if box_atoms is not None:
        options[ "box" ] = box_around ( coords[ list ( box_atoms[0] ) ], box_atoms[1] )
    result = detect_cavities ( coords, radii, **options )
    volume, ids, best = match_pocket ( result, reference_atoms, min_overlap )
    return { "volume": volume, "cavities": ids, "overlap": best, "n_cavities": len ( result[ "cavities" ] ),
             "total_volume": float ( sum ( c[ "volume" ] for c in result[ "cavities" ] ) ) }


def pocket_volume_trajectory ( frames, radii, reference_atoms, options = None, min_overlap = 0.3, processes = None ):
    """ Volume of one pocket (given by its lining atoms in a reference frame)
        in every frame: frames is a list of (n,3) coordinate arrays. One
        worker process per frame. Returns a list of dicts (volume, ids of
        the matched cavities, best overlap, all cavities' volume). """
    options = dict ( options or { } )
    options[ "lining" ] = True
    jobs = [ ( np.asarray ( f, dtype = np.float64 ), radii, options, list ( reference_atoms ), min_overlap ) for f in frames ]
    if len ( jobs ) <= 1 or processes == 1:
        return [ _pocket_frame_job ( j ) for j in jobs ]
    import multiprocessing
    with multiprocessing.Pool ( processes = min ( processes or multiprocessing.cpu_count ( ), len ( jobs ) ) ) as pool:
        return pool.map ( _pocket_frame_job, jobs )
