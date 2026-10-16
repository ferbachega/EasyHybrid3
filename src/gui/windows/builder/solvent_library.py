#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  solvent_library.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  [EN] 2026-10-03 user request: other solvents besides water, easy to add
#  ("talvez na forma de arquivos mol2").
#
#  The solvent library is the folder  src/gui/windows/builder/solvents/ :
#  ONE .mol2 FILE PER SOLVENT -- dropping a file there adds the solvent.
#  Two kinds of file are accepted:
#
#    * a PRE-EQUILIBRATED PERIODIC BOX: many copies of the molecule (one
#      SUBSTRUCTURE id per molecule, column 7 of the ATOM lines) and a
#      @<TRIPOS>CRYSIN section with the (orthorhombic) cell. Used directly as
#      the solvent template -- best quality (e.g. water_tip3p.mol2, or a box
#      exported from an MD run).
#    * a SINGLE MOLECULE: EasyHybrid builds the periodic template itself at
#      the solvent's density (lattice + random orientations, then random
#      insertion to reach the target count; see build_template_box()). The
#      result is cached in solvents/.cache/. Not equilibrated: minimise and
#      equilibrate (NPT) the solvated system.
#
#  Metadata: comment lines ("# key: value") at the top of the file --
#      # name: Methanol        label shown in the Solvate window
#      # residue: MOH          residue name of every solvent molecule
#      # density: 0.792        g/cm3 (required for single-molecule files)
#  Missing name -> the MOLECULE name line; missing residue -> the ATOM lines'
#  substructure name (else "SOL").
#
#  Atom names come from the file (made unique per molecule if needed), bond
#  orders from its BOND section ("ar" bonds Kekulized with the app's own
#  Wang & Case perceiver, as in structure_library.py).
# ============================================================================
import os
import hashlib
import numpy as np

SOLVENT_DIR = os.path.join ( os.path.dirname ( os.path.abspath ( __file__ ) ), "solvents" )
CACHE_DIR   = os.path.join ( SOLVENT_DIR, ".cache" )
AVOGADRO    = 6.02214076e23

_MOL2_BOND_ORDERS = { "1": 1, "2": 2, "3": 3, "ar": 1, "am": 1, "du": 1, "un": 1 }
_ATOMIC_MASS = { "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "P": 30.974,
                 "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.904, "Si": 28.085, "B": 10.81 }


class SolventError ( Exception ):
    pass


class Solvent:
    """ One solvent: the molecule (elements, names, residue, bonds) and its
        periodic template (molecules (N, n_atoms, 3), whole, centroids in
        [0, edges)). """
    def __init__ ( self ):
        self.key      = None        # file stem
        self.path     = None
        self.name     = None
        self.residue  = "SOL"
        self.density  = None        # g/cm3
        self.kind     = None        # "box" | "molecule"
        self.elements = [ ]
        self.names    = [ ]
        self.bonds    = [ ]         # [ (i, j, order) ] within ONE molecule
        self.template = None        # (N, n_atoms, 3)
        self.edges    = None        # (3,) orthorhombic template cell
        self.contact_any   = 1.4    # closest intermolecular contacts allowed at seams
        self.contact_heavy = 2.3
        self.notes    = [ ]

    @property
    def n_atoms ( self ):
        return len ( self.elements )

    @property
    def molar_mass ( self ):
        return float ( sum ( _ATOMIC_MASS.get ( e, 12.0 ) for e in self.elements ) )

    @property
    def molarity ( self ):
        """ mol/L of the pure solvent (55.5 for water). """
        density = self.density if self.density else self.template_density
        return density * 1000.0 / self.molar_mass

    @property
    def template_density ( self ):
        if self.template is None or self.edges is None:
            return None
        return len ( self.template ) * self.molar_mass / ( AVOGADRO * float ( np.prod ( self.edges ) ) * 1e-24 )

    @property
    def is_water ( self ):
        return sorted ( self.elements ) == [ "H", "H", "O" ]


# ----------------------------------------------------------------------------
#  mol2 reading
# ----------------------------------------------------------------------------
def _element ( atom_name, atom_type ):
    """ Element from the Tripos atom type ("C.3", "Cl", "O.co2"...), falling
        back to the atom name's letters. """
    base = atom_type.split ( "." )[0]
    for candidate in ( base, "".join ( c for c in atom_name if c.isalpha ( ) ) ):
        if not candidate:
            continue
        for length in ( 2, 1 ):
            symbol = candidate[:length].capitalize ( )
            if symbol in _ATOMIC_MASS or symbol in ( "Na", "K", "Li", "Mg", "Ca", "Zn", "Fe" ):
                return symbol
    return None


def read_header ( path ):
    """ {key: value} from the leading "# key: value" comment lines. """
    meta = { }
    with open ( path ) as handle:
        for line in handle:
            stripped = line.strip ( )
            if stripped.startswith ( "@<TRIPOS>" ):
                break
            if stripped.startswith ( "#" ) and ":" in stripped:
                key, value = stripped[1:].split ( ":", 1 )
                meta.setdefault ( key.strip ( ).lower ( ), value.strip ( ) )
    return meta


def _parse_mol2 ( path ):
    """ (meta, molecule_name, atoms, bonds, crysin)
        atoms: [ (name, element, x, y, z, subst_id, subst_name) ]
        bonds: [ (i, j, code) ] 0-based; crysin: (a, b, c, al, be, ga) or None """
    meta = read_header ( path )
    with open ( path ) as handle:
        lines = [ line.rstrip ( "\n" ) for line in handle ]
    atoms, bonds, crysin, molecule_name = [ ], [ ], None, None
    section = None
    n_atoms = n_bonds = None
    k = 0
    while k < len ( lines ):
        line = lines[k].strip ( )
        if line.startswith ( "@<TRIPOS>" ):
            section = line
            if section == "@<TRIPOS>MOLECULE":
                molecule_name = lines[k + 1].strip ( )
                counts = lines[k + 2].split ( )
                n_atoms = int ( counts[0] )
                n_bonds = int ( counts[1] ) if len ( counts ) > 1 else 0
                k += 3
                continue
        elif line and not line.startswith ( "#" ):
            fields = line.split ( )
            if section == "@<TRIPOS>ATOM" and len ( fields ) >= 6:
                element = _element ( fields[1], fields[5] )
                if element is None:
                    raise SolventError ( "{}: atom {!r} ({}) is not a recognised element.".format (
                            os.path.basename ( path ), fields[1], fields[5] ) )
                subst_id   = int ( fields[6] ) if len ( fields ) > 6 else 1
                subst_name = fields[7] if len ( fields ) > 7 else ""
                atoms.append ( ( fields[1], element, float ( fields[2] ), float ( fields[3] ), float ( fields[4] ),
                                 subst_id, subst_name ) )
            elif section == "@<TRIPOS>BOND" and len ( fields ) >= 4:
                bonds.append ( ( int ( fields[1] ) - 1, int ( fields[2] ) - 1, fields[3].lower ( ) ) )
            elif section == "@<TRIPOS>CRYSIN" and len ( fields ) >= 6 and crysin is None:
                crysin = tuple ( float ( v ) for v in fields[:6] )
        k += 1
    if n_atoms is None or len ( atoms ) != n_atoms:
        raise SolventError ( "{}: malformed mol2 (atom count).".format ( os.path.basename ( path ) ) )
    if n_bonds is not None and len ( bonds ) != n_bonds:
        raise SolventError ( "{}: malformed mol2 (bond count).".format ( os.path.basename ( path ) ) )
    return meta, molecule_name, atoms, bonds, crysin


def _resolve_orders ( elements, bonds ):
    """ [(i, j, order)]: explicit 1/2/3 kept, "ar" Kekulized (Wang & Case). """
    result = [ [ i, j, _MOL2_BOND_ORDERS.get ( code, 1 ), code == "ar" ] for ( i, j, code ) in bonds ]
    if any ( entry[3] for entry in result ):
        try:
            from vismol.core.bond_order_perception import perceive_bond_orders
            order_map, _ = perceive_bond_orders ( list ( elements ), [ ( i, j ) for ( i, j, _o, _a ) in result ] )
            for entry in result:
                if entry[3]:
                    key = ( min ( entry[0], entry[1] ), max ( entry[0], entry[1] ) )
                    entry[2] = int ( order_map.get ( key, 1 ) )
        except Exception:
            pass
    return [ ( i, j, order ) for ( i, j, order, _a ) in result ]


def _unique_names ( names, elements ):
    if len ( set ( names ) ) == len ( names ):
        return list ( names )
    counter, unique = { }, [ ]
    for element in elements:
        counter[element] = counter.get ( element, 0 ) + 1
        unique.append ( "{}{}".format ( element.upper ( ), counter[element] ) )
    return unique


def _make_whole_and_wrap ( molecules, edges ):
    """ Each molecule made whole (minimum image w.r.t. its first atom) and
        shifted so its centroid lies in [0, edges). """
    first = molecules[:, :1, :]
    delta = molecules - first
    delta -= edges * np.round ( delta / edges )
    whole = first + delta
    centroids = whole.mean ( axis = 1, keepdims = True )
    return whole - edges * np.floor ( centroids / edges )


# ----------------------------------------------------------------------------
#  library
# ----------------------------------------------------------------------------
def discover_solvents ( folder = None ):
    """ [ (key, label, path, kind) ] for every readable .mol2 in the
        library folder, water first, then alphabetically. Unreadable files
        are skipped (reported in the label list as nothing). """
    folder = folder or SOLVENT_DIR
    entries = [ ]
    if not os.path.isdir ( folder ):
        return entries
    for filename in sorted ( os.listdir ( folder ) ):
        if not filename.lower ( ).endswith ( ".mol2" ):
            continue
        path = os.path.join ( folder, filename )
        try:
            meta = read_header ( path )
            with open ( path ) as handle:
                is_box = "@<TRIPOS>CRYSIN" in handle.read ( )
        except Exception:
            continue
        key   = os.path.splitext ( filename )[0]
        label = meta.get ( "name" ) or key
        entries.append ( ( key, label, path, "box" if is_box else "molecule" ) )
    entries.sort ( key = lambda e: ( not e[0].startswith ( "water" ), e[1].lower ( ) ) )
    return entries


_loaded = { }


def load_solvent ( path, progress = None ):
    """ Solvent for a library file (cached in memory, and on disk for
        generated templates). Raises SolventError with a readable message. """
    stamp = ( os.path.getmtime ( path ), os.path.getsize ( path ) )
    cached = _loaded.get ( path )
    if cached is not None and cached[0] == stamp:
        return cached[1]

    meta, molecule_name, atoms, bonds, crysin = _parse_mol2 ( path )
    if not atoms:
        raise SolventError ( "{}: no atoms.".format ( os.path.basename ( path ) ) )
    solvent = Solvent ( )
    solvent.path = path
    solvent.key  = os.path.splitext ( os.path.basename ( path ) )[0]
    solvent.name = meta.get ( "name" ) or molecule_name or solvent.key
    if meta.get ( "density" ):
        try:
            solvent.density = float ( meta["density"].split ( )[0] )
        except ValueError:
            raise SolventError ( "{}: invalid density {!r}.".format ( os.path.basename ( path ), meta["density"] ) )

    subst_ids = [ a[5] for a in atoms ]
    groups = { }
    for index, sid in enumerate ( subst_ids ):
        groups.setdefault ( sid, [ ] ).append ( index )
    molecules_idx = list ( groups.values ( ) )
    first = molecules_idx[0]
    solvent.elements = [ atoms[i][1] for i in first ]
    solvent.names    = _unique_names ( [ atoms[i][0] for i in first ], solvent.elements )
    solvent.residue  = ( meta.get ( "residue" ) or atoms[first[0]][6] or "SOL" )[:4]
    position = { atom: k for k, atom in enumerate ( first ) }
    solvent.bonds = _resolve_orders ( solvent.elements,
                        [ ( position[i], position[j], code ) for ( i, j, code ) in bonds
                          if i in position and j in position ] )

    xyz = np.array ( [ a[2:5] for a in atoms ], dtype = np.float64 )
    if crysin is not None:
        solvent.kind = "box"
        a, b, c, alpha, beta, gamma = crysin
        if any ( abs ( angle - 90.0 ) > 1e-3 for angle in ( alpha, beta, gamma ) ):
            raise SolventError ( "{}: only orthorhombic boxes (90 degree angles) are supported.".format (
                    os.path.basename ( path ) ) )
        n = solvent.n_atoms
        if any ( len ( g ) != n for g in molecules_idx ) or \
           any ( [ atoms[i][1] for i in g ] != solvent.elements for g in molecules_idx ):
            raise SolventError ( "{}: every molecule of the box must have the same atoms in the same order "
                                 "(one SUBSTRUCTURE id per molecule).".format ( os.path.basename ( path ) ) )
        edges = np.array ( [ a, b, c ], dtype = np.float64 )
        molecules = np.stack ( [ xyz[g] for g in molecules_idx ] )
        solvent.template = _make_whole_and_wrap ( molecules, edges )
        solvent.edges = edges
        if solvent.density is None:
            solvent.density = solvent.template_density
    else:
        solvent.kind = "molecule"
        if len ( molecules_idx ) != 1:
            raise SolventError ( "{}: several molecules but no CRYSIN cell: give either ONE molecule or a "
                                 "periodic box with a @<TRIPOS>CRYSIN section.".format ( os.path.basename ( path ) ) )
        if not solvent.density:
            raise SolventError ( "{}: a single-molecule solvent needs its density -- add a line "
                                 "'# density: <g/cm3>' at the top of the file.".format ( os.path.basename ( path ) ) )
        solvent.template, solvent.edges = _cached_template ( solvent, xyz, progress )
    solvent.contact_any, solvent.contact_heavy = _template_contacts ( solvent )
    _loaded[path] = ( stamp, solvent )
    return solvent


def _template_contacts ( solvent ):
    """ Closest intermolecular contacts inside the template (minimum image),
        x 0.95 -- what counts as a clash at the seams of a tiled box (for the
        water template: 1.38 A any atom, 2.28 A O-O). """
    molecules = solvent.template
    edges     = solvent.edges
    flat      = molecules.reshape ( -1, 3 )
    owner     = np.repeat ( np.arange ( len ( molecules ) ), solvent.n_atoms )
    heavy     = np.tile ( np.array ( [ e != "H" for e in solvent.elements ] ), len ( molecules ) )
    best_any, best_heavy = 3.0, 4.0
    sample = np.arange ( 0, len ( molecules ), max ( 1, len ( molecules ) // 120 ) )
    for m in sample:
        atoms_m = molecules[m]
        delta = flat[None, :, :] - atoms_m[:, None, :]
        delta -= edges * np.round ( delta / edges )
        d = np.sqrt ( ( delta ** 2 ).sum ( axis = 2 ) )
        d[:, owner == m] = 99.0
        best_any = min ( best_any, float ( d.min ( ) ) )
        heavy_m = np.array ( [ e != "H" for e in solvent.elements ] )
        if heavy_m.any ( ) and heavy.any ( ):
            best_heavy = min ( best_heavy, float ( d[heavy_m][:, heavy].min ( ) ) )
    return max ( 1.2, 0.95 * best_any ), max ( 2.0, 0.95 * best_heavy )


# ----------------------------------------------------------------------------
#  template generation for single-molecule files
# ----------------------------------------------------------------------------
def _cached_template ( solvent, xyz, progress ):
    os.makedirs ( CACHE_DIR, exist_ok = True )
    digest = hashlib.sha1 ( open ( solvent.path, "rb" ).read ( ) ).hexdigest ( )[:12]
    cache = os.path.join ( CACHE_DIR, "{}_{}.npz".format ( solvent.key, digest ) )
    if os.path.isfile ( cache ):
        try:
            data = np.load ( cache )
            return data["template"], data["edges"]
        except Exception:
            pass
    template, edges = build_template_box ( xyz, solvent.elements, solvent.density, solvent.molar_mass,
                                           progress = progress )
    try:
        for stale in os.listdir ( CACHE_DIR ):            # older versions of this file
            if stale.startswith ( solvent.key + "_" ) and stale.endswith ( ".npz" ):
                os.remove ( os.path.join ( CACHE_DIR, stale ) )
        np.savez ( cache, template = template, edges = edges )
    except OSError:
        pass
    return template, edges


def _random_rotation ( rng ):
    q = rng.normal ( size = 4 )
    q /= np.linalg.norm ( q )
    w, x, y, z = q
    return np.array ( [ [ 1 - 2 * ( y * y + z * z ), 2 * ( x * y - z * w ), 2 * ( x * z + y * w ) ],
                        [ 2 * ( x * y + z * w ), 1 - 2 * ( x * x + z * z ), 2 * ( y * z - x * w ) ],
                        [ 2 * ( x * z - y * w ), 2 * ( y * z + x * w ), 1 - 2 * ( x * x + y * y ) ] ] )


def _rotation_about ( axis, angle ):
    axis = axis / np.linalg.norm ( axis )
    x, y, z = axis
    c, s_, t = np.cos ( angle ), np.sin ( angle ), 1.0 - np.cos ( angle )
    return np.array ( [ [ t * x * x + c,     t * x * y - s_ * z, t * x * z + s_ * y ],
                        [ t * x * y + s_ * z, t * y * y + c,     t * y * z - s_ * x ],
                        [ t * x * z - s_ * y, t * y * z + s_ * x, t * z * z + c     ] ] )


def build_template_box ( molecule_xyz, elements, density, molar_mass, edge = None, seed = 7,
                         d_any = 1.8, d_heavy = 2.8, max_sweeps = 600, progress = None ):
    """ Periodic cubic box of copies of one molecule at `density` (g/cm3),
        Packmol-style: ALL the molecules the density asks for are put on a
        lattice with random orientations, then overlaps are removed by rigid
        Monte Carlo moves (small translation + rotation of an overlapping
        molecule, accepted only when its overlap decreases; an occasional
        full re-orientation unsticks jammed ones). Overlap = sum over
        intermolecular atom pairs of max(0, limit - r)^2, limit = d_heavy for
        two heavy atoms, d_any otherwise (minimum image). Molecules still
        overlapping after `max_sweeps` are removed. Returns (template
        (N, n, 3), edges (3,)). """
    rng = np.random.default_rng ( seed )
    xyz = np.asarray ( molecule_xyz, dtype = np.float64 )
    xyz = xyz - xyz.mean ( axis = 0 )
    n_at = len ( xyz )
    radius = float ( np.sqrt ( ( xyz ** 2 ).sum ( axis = 1 ) ).max ( ) )
    if edge is None:
        edge = max ( 30.0, 6.0 * radius + 10.0 )
    target = max ( 1, int ( round ( density * 1e-24 * edge ** 3 * AVOGADRO / molar_mass ) ) )
    heavy  = np.array ( [ e != "H" for e in elements ] )

    per_side = int ( np.ceil ( target ** ( 1.0 / 3.0 ) ) )
    spacing  = edge / per_side
    sites = np.array ( [ ( i + 0.5, j + 0.5, k + 0.5 ) for i in range ( per_side )
                         for j in range ( per_side ) for k in range ( per_side ) ] ) * spacing
    sites = sites[rng.permutation ( len ( sites ) )[:target]]
    rotations = np.stack ( [ _random_rotation ( rng ) for _ in range ( target ) ] )
    centres   = sites.copy ( )
    atoms     = np.einsum ( "nij,aj->nai", rotations, xyz ) + centres[:, None, :]      # (N, n_at, 3)

    pair_limit = np.where ( heavy[:, None] & heavy[None, :], d_heavy, d_any )      # (n_at, n_at)
    cutoff_centres = 2.0 * radius + d_heavy

    def overlap ( m, coords ):
        """ Overlap of molecule m (at `coords`) with every OTHER molecule
            whose centre is close enough to touch it. """
        centre = coords.mean ( axis = 0 )
        d = centres - centre
        d -= edge * np.round ( d / edge )
        near = np.flatnonzero ( ( d * d ).sum ( axis = 1 ) < cutoff_centres ** 2 )
        near = near[near != m]
        if near.size == 0:
            return 0.0
        delta = atoms[near][:, None, :, :] - coords[None, :, None, :]           # (k, n_at, n_at, 3)
        delta -= edge * np.round ( delta / edge )
        r = np.sqrt ( ( delta ** 2 ).sum ( axis = 3 ) )
        excess = np.clip ( pair_limit[None, :, :] - r, 0.0, None )
        return float ( ( excess ** 2 ).sum ( ) )

    penalties = np.array ( [ overlap ( m, atoms[m] ) for m in range ( target ) ] )
    step_t, step_r = 0.35 * spacing / 3.0, np.radians ( 25.0 )
    for sweep in range ( max_sweeps ):
        bad = np.flatnonzero ( penalties > 1e-9 )
        if bad.size == 0:
            break
        for m in rng.permutation ( bad ):
            if rng.random ( ) < 0.05:                     # full re-orientation
                new_rot = _random_rotation ( rng )
                new_centre = centres[m] + rng.normal ( 0.0, step_t, 3 )
            else:
                new_rot = _rotation_about ( rng.normal ( size = 3 ), rng.uniform ( -step_r, step_r ) ) @ rotations[m]
                new_centre = centres[m] + rng.normal ( 0.0, step_t, 3 )
            new_atoms = xyz @ new_rot.T + new_centre
            new_pen = overlap ( m, new_atoms )
            if new_pen < penalties[m]:
                # the neighbours' penalties change too: recompute those near m
                atoms[m], centres[m], rotations[m] = new_atoms, new_centre, new_rot
                d = centres - new_centre
                d -= edge * np.round ( d / edge )
                near = np.flatnonzero ( np.sqrt ( ( d ** 2 ).sum ( axis = 1 ) ) < 2.0 * radius + d_heavy + step_t * 4 )
                for k in near:
                    penalties[k] = overlap ( k, atoms[k] )
        if progress is not None and sweep % 10 == 0:
            progress ( 1.0 - bad.size / target )
    # anything still overlapping: drop the worst offender until none is left
    keep = np.ones ( target, dtype = bool )
    if ( penalties > 1e-9 ).any ( ):
        while True:
            remaining = np.flatnonzero ( keep )
            current = { k: _overlap_subset ( k, atoms, remaining, heavy, d_any, d_heavy, edge )
                        for k in remaining if penalties[k] > 1e-9 }
            bad = { k: v for k, v in current.items ( ) if v > 1e-9 }
            if not bad:
                break
            keep[max ( bad, key = bad.get )] = False
    template = atoms[keep]
    edges = np.full ( 3, edge )
    return _make_whole_and_wrap ( template, edges ), edges


def _overlap_subset ( m, atoms, subset, heavy, d_any, d_heavy, edge ):
    others = atoms[[ k for k in subset if k != m ]].reshape ( -1, 3 )
    if len ( others ) == 0:
        return 0.0
    heavy_o = np.tile ( heavy, len ( others ) // len ( heavy ) )
    delta = others[None, :, :] - atoms[m][:, None, :]
    delta -= edge * np.round ( delta / edge )
    r = np.sqrt ( ( delta ** 2 ).sum ( axis = 2 ) )
    limit = np.where ( heavy[:, None] & heavy_o[None, :], d_heavy, d_any )
    return float ( ( np.clip ( limit - r, 0.0, None ) ** 2 ).sum ( ) )


def quick_info ( path ):
    """ Cheap properties for the Solvate window's estimate, WITHOUT building
        a template: dict(name, residue, kind, density, molar_mass, n_atoms,
        number_density (molecules/A3), cached (template already available)).
        Raises SolventError for an unreadable file. """
    meta, molecule_name, atoms, bonds, crysin = _parse_mol2 ( path )
    if not atoms:
        raise SolventError ( "{}: no atoms.".format ( os.path.basename ( path ) ) )
    first_id = atoms[0][5]
    first    = [ a for a in atoms if a[5] == first_id ]
    mass     = float ( sum ( _ATOMIC_MASS.get ( a[1], 12.0 ) for a in first ) )
    density  = None
    if meta.get ( "density" ):
        try:
            density = float ( meta["density"].split ( )[0] )
        except ValueError:
            density = None
    n_molecules = len ( set ( a[5] for a in atoms ) )
    if density is None and crysin is not None:
        volume  = crysin[0] * crysin[1] * crysin[2]
        density = n_molecules * mass / ( AVOGADRO * volume * 1e-24 )
    key = os.path.splitext ( os.path.basename ( path ) )[0]
    cached = path in _loaded
    if not cached and crysin is None and os.path.isdir ( CACHE_DIR ):
        cached = any ( name.startswith ( key + "_" ) for name in os.listdir ( CACHE_DIR ) )
    return dict ( name = meta.get ( "name" ) or molecule_name or key,
                  residue = ( meta.get ( "residue" ) or first[0][6] or "SOL" )[:4],
                  kind = "box" if crysin is not None else "molecule",
                  density = density, molar_mass = mass, n_atoms = len ( first ),
                  number_density = ( density * AVOGADRO * 1e-24 / mass ) if density else None,
                  cached = crysin is not None or cached )


def add_solvent_file ( source, name = None, residue = None, density = None, folder = None ):
    """ Copies a .mol2 file into the library, writing/overriding the "# name",
        "# residue" and "# density" header lines with the values given.
        Returns the new path (a numbered suffix avoids overwriting). """
    folder = folder or SOLVENT_DIR
    os.makedirs ( folder, exist_ok = True )
    _parse_mol2 ( source )                                     # readable?
    with open ( source ) as handle:
        text = handle.read ( )
    lines = text.split ( "\n" )
    start = next ( ( k for k, line in enumerate ( lines ) if line.strip ( ).startswith ( "@<TRIPOS>" ) ), 0 )
    header, body = lines[:start], lines[start:]
    overrides = { "name": name, "residue": residue,
                  "density": ( "{:.4f}".format ( density ) if density else None ) }
    kept = [ line for line in header
             if not ( line.strip ( ).startswith ( "#" ) and ":" in line
                      and line.strip ( )[1:].split ( ":", 1 )[0].strip ( ).lower ( ) in
                      [ k for k, v in overrides.items ( ) if v ] ) ]
    new_header = [ "# {}: {}".format ( k, v ) for k, v in overrides.items ( ) if v ] + kept
    stem = os.path.splitext ( os.path.basename ( source ) )[0].replace ( " ", "_" )
    target = os.path.join ( folder, stem + ".mol2" )
    n = 2
    while os.path.exists ( target ):
        target = os.path.join ( folder, "{}_{}.mol2".format ( stem, n ) )
        n += 1
    with open ( target, "w" ) as handle:
        handle.write ( "\n".join ( new_header + body ) )
    return target
