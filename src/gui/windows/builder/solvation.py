#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  solvation.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  [EN] 2026-10-03 user request: a solvation tool written from scratch (the
#  pDynamo route -- AddCounterIons + SolvateSystemBySuperposition -- "has not
#  worked well"), as a Builder tool, independent of any force field.
#
#  Pure numpy, no GUI and no pDynamo: takes solute coordinates + total charge,
#  returns where the solute must be moved to, the cell, the waters and the
#  ions. builder_solvation.py turns that into Builder atoms/residues.
#
#  Method (the tleap / VMD-solvate / gmx-solvate approach):
#    1. A pre-equilibrated periodic water box (solvents/water_tip3p_box.xyz:
#       734 rigid TIP3P waters, 28 A cube, 1.00 g/cm3) is replicated until it
#       covers the requested region. Inside the region the copies join
#       seamlessly (the template is itself periodic).
#    2. Waters whose oxygen falls outside the region are dropped; for a
#       periodic cell, waters that clash with the periodic image of a water on
#       the opposite face (the region edge is generally not a multiple of the
#       template edge) are dropped too.
#    3. Waters with any atom closer than `closeness` to any solute atom
#       (periodic images included) are dropped.
#    4. Ions: counter-ions to neutralise the solute's formal charge, plus
#       salt pairs for the requested concentration (n = c * N_water / 55.5,
#       the usual "per water" estimate). Each ion REPLACES a randomly chosen
#       water whose oxygen is at least `ion_solute_distance` from the solute
#       and `ion_ion_distance` from every ion already placed (minimum image);
#       both limits are relaxed step by step when the box is too crowded.
#
#  Shapes (region_for_solute()):
#    "cube"    periodic cube, edge = largest solute extent + 2 x margin
#    "box"     periodic rectangular box, each edge = solute extent along that
#              axis + 2 x margin (optionally after aligning the solute's
#              principal axes with x, y, z, which shrinks boxes around
#              elongated molecules)
#    "sphere"  non-periodic droplet (QM/MM water spheres, cluster models):
#              radius = farthest solute atom from the centre + margin. No
#              cell and no seams (the tiled template is cut by the sphere);
#              ions are also kept away from the droplet surface.
#  The truncated octahedron is planned (triclinic cell + its own minimum
#  image).
# ============================================================================
import os
import numpy as np

# monovalent ions offered by the tool: label -> (element, residue/atom name, charge)
CATIONS = { "Na+": ( "Na", "NA", +1 ), "K+": ( "K", "K", +1 ), "Li+": ( "Li", "LI", +1 ) }
ANIONS  = { "Cl-": ( "Cl", "CL", -1 ), "Br-": ( "Br", "BR", -1 ) }

SHAPES_IMPLEMENTED = ( "cube", "box", "sphere" )
WATER_MOLARITY     = 55.5          # mol/L of pure water


class SolvationError ( Exception ):
    pass


DEFAULT_SOLVENT_FILE = "water_tip3p.mol2"


def default_solvent ( ):
    """ The library's water (solvents/water_tip3p.mol2). """
    try:
        from gui.windows.builder.solvent_library import load_solvent, SOLVENT_DIR
    except ImportError:                          # used as a plain module (no gui package)
        from solvent_library import load_solvent, SOLVENT_DIR
    return load_solvent ( os.path.join ( SOLVENT_DIR, DEFAULT_SOLVENT_FILE ) )


def load_water_template ( ):
    """ (waters (N, 3, 3), edge) -- kept for older callers/tests. """
    water = default_solvent ( )
    return water.template, float ( water.edges[0] )


# ----------------------------------------------------------------------------
#  neighbour search (cell list, exact distances)
# ----------------------------------------------------------------------------
def _close_mask ( points, reference, cutoff ):
    """ Boolean mask: True where a point has a reference point closer than
        `cutoff`. Cell list on a `cutoff`-sized grid, exact distances. """
    points    = np.asarray ( points,    dtype = np.float64 ).reshape ( -1, 3 )
    reference = np.asarray ( reference, dtype = np.float64 ).reshape ( -1, 3 )
    mask = np.zeros ( len ( points ), dtype = bool )
    if len ( points ) == 0 or len ( reference ) == 0:
        return mask
    origin  = np.minimum ( points.min ( axis = 0 ), reference.min ( axis = 0 ) )
    ref_key = np.floor ( ( reference - origin ) / cutoff ).astype ( np.int64 )
    pt_key  = np.floor ( ( points    - origin ) / cutoff ).astype ( np.int64 )
    cells = { }
    order = np.lexsort ( ref_key.T[::-1] )
    sorted_keys = ref_key[order]
    boundaries = np.flatnonzero ( np.any ( np.diff ( sorted_keys, axis = 0 ) != 0, axis = 1 ) ) + 1
    for chunk in np.split ( order, boundaries ):
        cells[tuple ( ref_key[chunk[0]] )] = chunk
    cutoff2 = cutoff * cutoff
    unique_keys, inverse = np.unique ( pt_key, axis = 0, return_inverse = True )
    inverse = inverse.reshape ( -1 )
    offsets = [ ( i, j, k ) for i in ( -1, 0, 1 ) for j in ( -1, 0, 1 ) for k in ( -1, 0, 1 ) ]
    point_groups = np.split ( np.argsort ( inverse, kind = "stable" ),
                              np.cumsum ( np.bincount ( inverse ) )[:-1] )
    for key_index, key in enumerate ( unique_keys ):
        kx, ky, kz = key
        neighbours = [ cells[c] for c in ( ( kx + i, ky + j, kz + k ) for ( i, j, k ) in offsets ) if c in cells ]
        if not neighbours:
            continue
        ref = reference[np.concatenate ( neighbours )]
        idx = point_groups[key_index]
        d2  = ( ( points[idx, None, :] - ref[None, :, :] ) ** 2 ).sum ( axis = 2 )
        mask[idx] = ( d2 < cutoff2 ).any ( axis = 1 )
    return mask


def _periodic_images_near_faces ( xyz, cell_lengths, margin ):
    """ Copies of the atoms within `margin` of a face of the orthorhombic
        cell [0, L), shifted by +/-L, so a plain (non-periodic) neighbour
        search also sees the periodic neighbours across that face. """
    xyz = np.asarray ( xyz, dtype = np.float64 ).reshape ( -1, 3 )
    L   = np.asarray ( cell_lengths, dtype = np.float64 )
    images = [ ]
    for shift in ( ( i, j, k ) for i in ( -1, 0, 1 ) for j in ( -1, 0, 1 ) for k in ( -1, 0, 1 ) ):
        if shift == ( 0, 0, 0 ):
            continue
        s = np.array ( shift, dtype = np.float64 )
        moved = xyz + s * L
        keep = np.all ( ( moved > -margin ) & ( moved < L + margin ), axis = 1 )
        if keep.any ( ):
            images.append ( moved[keep] )
    return np.concatenate ( images ) if images else np.zeros ( ( 0, 3 ) )


def _minimum_image_distance ( a, b, cell_lengths ):
    d = np.asarray ( a ) - np.asarray ( b )
    if cell_lengths is not None:
        L = np.asarray ( cell_lengths )
        d = d - L * np.round ( d / L )
    return np.sqrt ( ( d * d ).sum ( axis = -1 ) )


# ----------------------------------------------------------------------------
#  regions
# ----------------------------------------------------------------------------
class Region:
    """ periodic: lengths (orthorhombic cell [0, L)) ; sphere: centre +
        radius. `solute_xyz` = solute coordinates after the move (and the
        optional rotation) the shape asks for. """
    def __init__ ( self, shape ):
        self.shape      = shape
        self.periodic   = shape in ( "cube", "box" )
        self.lengths    = None
        self.centre     = None
        self.radius     = None
        self.solute_xyz = None

    @property
    def volume ( self ):
        if self.periodic:
            return float ( np.prod ( self.lengths ) )
        return 4.0 / 3.0 * np.pi * self.radius ** 3

    def describe ( self ):
        if self.periodic:
            return "box {:.1f} x {:.1f} x {:.1f} \u00c5".format ( *self.lengths )
        return "sphere of radius {:.1f} \u00c5".format ( self.radius )


def principal_axes_rotation ( xyz ):
    """ Rotation matrix (rows = new x, y, z) aligning the largest spread of
        `xyz` with x, the next with y (unweighted PCA). """
    xyz = np.asarray ( xyz, dtype = np.float64 ).reshape ( -1, 3 )
    if len ( xyz ) < 3:
        return np.eye ( 3 )
    centred = xyz - xyz.mean ( axis = 0 )
    values, vectors = np.linalg.eigh ( centred.T @ centred )
    axes = vectors[:, np.argsort ( values )[::-1]].T
    if np.linalg.det ( axes ) < 0:                 # keep it a proper rotation
        axes[2] *= -1.0
    return axes


def region_for_solute ( solute_xyz, shape = "cube", padding = 10.0, box_size = None, align = False ):
    """ Region for a shape (see the module docstring). box_size: fixed size
        instead of the margin -- edge (cube), (a, b, c) (box) or radius
        (sphere). Periodic regions span [0, L) with the solute's bounding-box
        centre moved to the middle (vismol draws cells from the origin); a
        sphere keeps the solute where it is. """
    if shape not in SHAPES_IMPLEMENTED:
        raise SolvationError ( "Shape '{}' is not implemented yet (available: {}).".format (
                shape, ", ".join ( SHAPES_IMPLEMENTED ) ) )
    region = Region ( shape )
    xyz = np.asarray ( solute_xyz, dtype = np.float64 ).reshape ( -1, 3 )
    if align and shape == "box" and len ( xyz ) >= 3:
        centre0 = xyz.mean ( axis = 0 )
        xyz = ( xyz - centre0 ) @ principal_axes_rotation ( xyz ).T + centre0
    if len ( xyz ):
        low, high = xyz.min ( axis = 0 ), xyz.max ( axis = 0 )
    else:
        low = high = np.zeros ( 3 )
    centre = 0.5 * ( low + high )

    if shape == "sphere":
        if box_size is not None and np.ndim ( box_size ) == 0 and box_size > 0:
            radius = float ( box_size )
        else:
            reach  = float ( np.sqrt ( ( ( xyz - centre ) ** 2 ).sum ( axis = 1 ) ).max ( ) ) if len ( xyz ) else 0.0
            radius = reach + padding
        if radius <= 3.0:
            raise SolvationError ( "The sphere radius ({:.2f} A) is too small.".format ( radius ) )
        region.centre, region.radius, region.solute_xyz = centre, radius, xyz
        return region

    if shape == "cube":
        if box_size is not None and np.ndim ( box_size ) == 0 and box_size > 0:
            lengths = np.full ( 3, float ( box_size ) )
        else:
            lengths = np.full ( 3, float ( ( high - low ).max ( ) + 2.0 * padding ) )
    else:   # box
        if box_size is not None and np.ndim ( box_size ) == 1 and len ( box_size ) == 3 and min ( box_size ) > 0:
            lengths = np.array ( box_size, dtype = np.float64 )
        else:
            lengths = ( high - low ) + 2.0 * padding
    if lengths.min ( ) <= 3.0:
        raise SolvationError ( "The box ({:.2f} x {:.2f} x {:.2f} A) is too small.".format ( *lengths ) )
    region.lengths    = lengths
    region.solute_xyz = xyz + ( 0.5 * lengths - centre )
    return region


# ----------------------------------------------------------------------------
#  waters
# ----------------------------------------------------------------------------
def _tile ( solvent, first, last ):
    """ Copies of the solvent template for the integer block range
        [first, last] (inclusive) in each axis. """
    template, edges = solvent.template, solvent.edges
    blocks = [ template + np.array ( [ i, j, k ] ) * edges
               for i in range ( first[0], last[0] + 1 )
               for j in range ( first[1], last[1] + 1 )
               for k in range ( first[2], last[2] + 1 ) ]
    return np.concatenate ( blocks ) if blocks else np.zeros ( ( 0, solvent.n_atoms, 3 ) )


def _fill_box ( lengths, solvent ):
    """ Solvent molecules (N, n_atoms, 3) filling the periodic orthorhombic
        cell [0, L): the template tiled, molecules kept by centroid. Where L
        is not a multiple of the template edge the two sides of each periodic
        seam are independent liquid surfaces, and removing their clashes
        costs ~4-6 % of the density (water: 0.94-0.96 g/cm3; tleap's
        solvateBox behaves the same) -- minimisation + NPT restores it. The
        clash thresholds are the closest contacts found inside the template
        itself (solvent.contact_any / contact_heavy). """
    n_copies  = np.ceil ( lengths / solvent.edges ).astype ( int )
    molecules = _tile ( solvent, ( 0, 0, 0 ), n_copies - 1 )
    centroids = molecules.mean ( axis = 1 )
    molecules = molecules[np.all ( ( centroids >= 0.0 ) & ( centroids < lengths ), axis = 1 )]
    if np.any ( np.abs ( lengths / solvent.edges - np.round ( lengths / solvent.edges ) ) > 1e-6 ):
        molecules = _remove_seam_clashes ( molecules, lengths, solvent )
    return molecules


def _fill_box_with_water ( lengths, seam_oo = None, seam_any = None ):
    """ Water-only shortcut (older callers/tests). """
    return _fill_box ( np.asarray ( lengths, dtype = np.float64 ), default_solvent ( ) )


def _fill_sphere ( centre, radius, solvent ):
    """ Molecules whose centroid lies inside the sphere: the periodic
        template tiled over the sphere's bounding cube (seamless anywhere,
        the template being periodic) and cut by the sphere. """
    first = np.floor ( ( centre - radius ) / solvent.edges ).astype ( int )
    last  = np.floor ( ( centre + radius ) / solvent.edges ).astype ( int )
    molecules = _tile ( solvent, first, last )
    centroids = molecules.mean ( axis = 1 )
    return molecules[( ( centroids - centre ) ** 2 ).sum ( axis = 1 ) <= radius * radius]


def _remove_seam_clashes ( molecules, lengths, solvent ):
    """ Drops molecules that, THROUGH a periodic boundary, come closer than
        the template's own closest contacts to a molecule kept before them.
        Pairs inside the box are never touched (the template already has good
        contacts there). """
    centroids = molecules.mean ( axis = 1 )
    reach  = float ( np.sqrt ( ( ( solvent.template - solvent.template.mean ( axis = 1, keepdims = True ) ) ** 2
                                 ).sum ( axis = 2 ) ).max ( ) )
    margin = 2.0 * reach + solvent.contact_heavy
    heavy  = np.array ( [ e != "H" for e in solvent.elements ] )
    limit  = np.where ( heavy[:, None] & heavy[None, :], solvent.contact_heavy, solvent.contact_any )
    near_face  = np.any ( ( centroids < margin ) | ( centroids > lengths - margin ), axis = 1 )
    candidates = np.flatnonzero ( near_face )
    removed = np.zeros ( len ( molecules ), dtype = bool )
    for position, a in enumerate ( candidates ):
        if removed[a]:
            continue
        others = candidates[position + 1:]
        others = others[~removed[others]]
        if others.size == 0:
            continue
        d_image  = _minimum_image_distance ( centroids[a], centroids[others], lengths )
        d_direct = np.sqrt ( ( ( centroids[a] - centroids[others] ) ** 2 ).sum ( axis = 1 ) )
        close = others[( d_image < d_direct - 1e-6 ) & ( d_image < margin )]
        for b in close:
            d = _minimum_image_distance ( molecules[a][:, None, :], molecules[b][None, :, :], lengths )
            if ( d < limit ).any ( ):
                removed[b] = True
    return molecules[~removed]


# ----------------------------------------------------------------------------
#  ions
# ----------------------------------------------------------------------------
def ion_counts ( solute_charge, n_waters, concentration = 0.15, neutralize = True, molarity = WATER_MOLARITY ):
    """ (n_cations, n_anions) for monovalent ions. concentration in mol/L;
        n_waters = number of solvent molecules, molarity = mol/L of the pure
        solvent (55.5 for water): n_pairs = c * N / molarity. """
    n_pairs = int ( round ( max ( 0.0, concentration ) * n_waters / molarity ) )
    q = int ( round ( solute_charge ) )
    n_cat = n_an = n_pairs
    if neutralize:
        if q > 0:
            n_an  += q
        elif q < 0:
            n_cat += -q
    return n_cat, n_an


def _place_ions ( waters, solute_xyz, lengths, n_ions, ion_solute_distance, ion_ion_distance, rng,
                  eligible = None ):
    """ Indices of the waters to replace, in placement order. lengths=None:
        non-periodic (no images). eligible: optional boolean mask of the
        waters that may be replaced at all (e.g. not on a droplet surface);
        dropped too if there are not enough of them. """
    if n_ions <= 0:
        return [ ]
    oxygens = waters.mean ( axis = 1 )          # molecule centroids (the O for water: ~0.07 A off)
    if n_ions > len ( waters ):
        raise SolvationError ( "Not enough waters ({}) to place {} ions.".format ( len ( waters ), n_ions ) )
    reference = np.asarray ( solute_xyz, dtype = np.float64 ).reshape ( -1, 3 )
    d_sol, d_ion = float ( ion_solute_distance ), float ( ion_ion_distance )
    chosen = [ ]
    while True:
        if len ( reference ):
            ref = reference if lengths is None else \
                  np.concatenate ( [ reference, _periodic_images_near_faces ( reference, lengths, d_sol ) ] )
            allowed = ~_close_mask ( oxygens, ref, d_sol ) if d_sol > 0 else np.ones ( len ( waters ), bool )
        else:
            allowed = np.ones ( len ( waters ), dtype = bool )
        if eligible is not None:
            allowed &= eligible
        allowed[chosen] = False
        candidates = np.flatnonzero ( allowed )
        rng.shuffle ( candidates )
        for c in candidates:
            if len ( chosen ) == n_ions:
                break
            if chosen and _minimum_image_distance ( oxygens[c], oxygens[chosen], lengths ).min ( ) < d_ion:
                continue
            chosen.append ( int ( c ) )
        if len ( chosen ) == n_ions:
            return chosen
        if d_sol <= 0.0 and d_ion <= 0.0:
            raise SolvationError ( "Could not place {} ions.".format ( n_ions ) )
        # too crowded: relax both limits a little and keep what is placed
        if eligible is not None and d_sol <= 0.0 and d_ion <= 0.0:
            eligible = None
            d_sol, d_ion = float ( ion_solute_distance ), float ( ion_ion_distance )
            continue
        d_sol = max ( 0.0, d_sol - 0.5 )
        d_ion = max ( 0.0, d_ion - 0.5 )


# ----------------------------------------------------------------------------
#  driver
# ----------------------------------------------------------------------------
class SolvationResult:
    def __init__ ( self ):
        self.translation = np.zeros ( 3 )     # add to the solute coordinates
        self.cell        = None               # (a, b, c, alpha, beta, gamma) or None
        self.sphere      = None               # (centre, radius) for a droplet
        self.region      = None
        self.solute_xyz  = None               # solute coordinates after the move
        self.solvent     = None               # solvent_library.Solvent
        self.waters      = np.zeros ( ( 0, 3, 3 ) )
        self.ions        = [ ]                # [ (label, element, name, charge, xyz) ]
        self.solute_charge = 0
        self.removed_by_solute = 0
        self.notes       = [ ]

    @property
    def n_waters ( self ):
        return len ( self.waters )

    @property
    def molecules ( self ):
        """ Solvent molecules (N, n_atoms, 3) -- `waters` is the same array
            (name kept from the water-only version). """
        return self.waters

    def summary ( self ):
        n_cat = sum ( 1 for ion in self.ions if ion[3] > 0 )
        n_an  = sum ( 1 for ion in self.ions if ion[3] < 0 )
        total = self.solute_charge + sum ( ion[3] for ion in self.ions )
        if self.solvent is None or self.solvent.is_water:
            what = "waters"
        else:
            what = "{} molecules".format ( self.solvent.name )
        text = "{} {}, {} cations, {} anions".format ( self.n_waters, what, n_cat, n_an )
        if self.cell is not None:
            text += "; cell {:.2f} x {:.2f} x {:.2f} A".format ( *self.cell[:3] )
        elif self.sphere is not None:
            text += "; non-periodic sphere, radius {:.2f} A".format ( self.sphere[1] )
        text += "; solute charge {:+d}, total charge {:+d}".format ( int ( self.solute_charge ), int ( total ) )
        return text


def solvate ( solute_xyz, solute_charge = 0, shape = "cube", padding = 10.0, box_size = None,
              closeness = 2.4, neutralize = True, concentration = 0.15, cation = "Na+", anion = "Cl-",
              ion_solute_distance = 5.0, ion_ion_distance = 5.0, align = False, ion_surface_distance = 4.0,
              seed = None, solvent = None ):
    """ See the module docstring. Coordinates in A, concentration in mol/L.
        Returns a SolvationResult; result.solute_xyz holds the solute's new
        coordinates (moved into the box / rotated when aligned; unchanged
        for a sphere). box_size: edge (cube), (a, b, c) (box) or radius
        (sphere) instead of the margin. ion_surface_distance (sphere only):
        ions stay at least this far inside the droplet surface. solvent: a
        solvent_library.Solvent (default: the library's TIP3P water). """
    if cation not in CATIONS:
        raise SolvationError ( "Unknown cation '{}'.".format ( cation ) )
    if anion not in ANIONS:
        raise SolvationError ( "Unknown anion '{}'.".format ( anion ) )
    solute_xyz = np.asarray ( solute_xyz, dtype = np.float64 ).reshape ( -1, 3 )
    if solvent is None:
        solvent = default_solvent ( )
    if solvent.template is None or len ( solvent.template ) == 0:
        raise SolvationError ( "The solvent '{}' has no template molecules.".format ( solvent.name ) )
    result = SolvationResult ( )
    result.solvent = solvent
    result.solute_charge = int ( round ( solute_charge ) )

    region = region_for_solute ( solute_xyz, shape, padding, box_size, align )
    result.region     = region
    result.solute_xyz = region.solute_xyz
    moved = region.solute_xyz
    if len ( moved ):
        result.translation = moved.mean ( axis = 0 ) - solute_xyz.mean ( axis = 0 )

    if region.periodic:
        lengths = region.lengths
        result.cell = ( float ( lengths[0] ), float ( lengths[1] ), float ( lengths[2] ), 90.0, 90.0, 90.0 )
        waters = _fill_box ( lengths, solvent )
        ref = np.concatenate ( [ moved, _periodic_images_near_faces ( moved, lengths, closeness ) ] ) if len ( moved ) else moved
    else:
        lengths = None
        result.sphere = ( region.centre.copy ( ), float ( region.radius ) )
        waters = _fill_sphere ( region.centre, region.radius, solvent )
        ref = moved
    if len ( moved ) and len ( waters ):
        clash = _close_mask ( waters.reshape ( -1, 3 ), ref, closeness ).reshape ( len ( waters ), -1 ).any ( axis = 1 )
        result.removed_by_solute = int ( clash.sum ( ) )
        waters = waters[~clash]

    eligible = None
    if not region.periodic and len ( waters ):
        reach = np.sqrt ( ( ( waters.mean ( axis = 1 ) - region.centre ) ** 2 ).sum ( axis = 1 ) )
        eligible = reach <= region.radius - ion_surface_distance

    rng = np.random.default_rng ( seed )
    n_cat, n_an = ion_counts ( result.solute_charge, len ( waters ), concentration, neutralize, solvent.molarity )
    # alternate cations and anions in the placement order so neither kind
    # gets only the leftovers of a crowded box
    order  = [ ]
    cats, ans = [ cation ] * n_cat, [ anion ] * n_an
    while cats or ans:
        if cats: order.append ( cats.pop ( ) )
        if ans:  order.append ( ans.pop ( ) )
    replaced = _place_ions ( waters, moved, lengths, len ( order ), ion_solute_distance, ion_ion_distance,
                             rng, eligible )
    for label, w in zip ( order, replaced ):
        element, name, charge = CATIONS.get ( label ) or ANIONS[label]
        result.ions.append ( ( label, element, name, charge, waters[w].mean ( axis = 0 ) ) )
    if replaced:
        keep = np.ones ( len ( waters ), dtype = bool )
        keep[replaced] = False
        waters = waters[keep]
    result.waters = waters
    return result
