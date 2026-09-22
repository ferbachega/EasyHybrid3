#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  structure_library.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Builder structure library -- user's own request: "adicionar estruturas
#  prontas (de uma biblioteca de arquivos mol2), ao clicar no background, o
#  easyhybrid adiciona a estrutura. semelhante a ferramenta de fragmentos,
#  mas nao tem necessidade de ter um H como referencia." A COMPLETE, ready-
#  made molecule (a ligand, a solvent molecule, a cofactor, whatever .mol2
#  file the user drops into the library folder) gets placed as a WHOLE onto
#  empty space in the 3D view -- see atom_ops.add_structure_at_position() /
#  click_mode.handle_click_to_add_structure().
#
#  [EN] Deliberately a SEPARATE module from fragment_library.py, not a
#  parameterised variant of it, because the two file conventions have
#  genuinely different semantics:
#    - fragment_library.load_fragment() REQUIRES exactly one atom left one
#      bond short of its standard valence (the attachment point) and
#      REJECTS aromatic ("ar") bonds outright (no Kekulization step exists
#      there to turn "ar" into a concrete order -- see that module's own
#      top-of-file note).
#    - load_structure() here expects a chemically COMPLETE molecule (every
#      atom already at its full valence -- nothing to check/enforce, unlike
#      a fragment's intentionally-open attachment atom) and DOES accept
#      aromatic bonds (real structure files -- ligands, drug-like molecules,
#      cofactors -- commonly have aromatic rings): an "ar" bond is recorded
#      as order=1 plus a separate aromatic flag, mirroring the Builder's own
#      "Aromatic" bond-order tool convention (see atom_ops.add_bond()'s own
#      docstring for why pDynamo3's own Bond model keeps order and
#      aromaticity independent).
#    - Element resolution here goes through vismol's own full PeriodicTable
#      (every element it knows, not just fragment_library.py's small,
#      valence-checked subset) -- a complete structure's bonds are trusted
#      as declared in the file, there's no per-atom valence validation to
#      restrict the element set for in the first place.
#
#  Reuses the exact same "walk base_dir's immediate subdirectories, one
#  category per subfolder, flat list of .mol2 files inside each" discovery
#  convention fragment_library.discover_fragment_files() already
#  established, just pointed at this module's own "structures" subfolder by
#  default, so a structure library is just as easy to grow by hand (drop a
#  .mol2 file into a category folder, no code change, no restart needed --
#  structure_library_window.py re-scans on every open).
# ============================================================================
import os
from vismol.utils.elements import PeriodicTable


# [EN] Every element symbol vismol's own PeriodicTable knows, keyed by its
# UPPER-CASE form -> proper-case symbol (e.g. "CL" -> "Cl") -- built once at
# import time. Unlike fragment_library.py's own small, hand-picked
# _SUPPORTED_ELEMENTS (restricted to what its valence-checking logic can
# reason about), a complete structure's bonds are trusted as declared in the
# file -- any element the rest of the app already knows how to render/color
# is accepted here.
_PERIODIC_TABLE_SYMBOLS = { symbol.upper ( ): symbol for symbol in PeriodicTable ( ).elements_by_symbol.keys ( ) }

# [EN] Tripos MOL2 bond-type codes this loader can represent. Unlike
# fragment_library.py's _MOL2_BOND_ORDERS, "ar" (aromatic) IS supported here
# -- see this module's own top-of-file note for why. "du"/"un"
# (dummy/unknown) and anything else are still rejected outright (nothing
# meaningful to fall back to for those).
_MOL2_BOND_ORDERS    = { '1': 1, '2': 2, '3': 3, 'am': 1, 'ar': 1 }
_MOL2_AROMATIC_CODES = { 'ar' }


class StructureError ( Exception ):
    """ Raised by load_structure() for any malformed/unsupported structure
    file -- always carries a short, specific, user-facing message, same
    convention as fragment_library.FragmentError. """
    pass


def discover_structure_files ( base_dir = None ):
    """ [EN] Same shape/contract as fragment_library.discover_fragment_
    files() -- walks base_dir's IMMEDIATE subdirectories (each one a
    category, e.g. "ligands"/"solvents"/"cofactors") and returns every
    ".mol2" file found directly inside each one.

    base_dir defaults to this module's own "structures" subdirectory
    (a sibling of fragment_library.py's own "fragments" folder) -- ships
    empty; the whole point of this feature is that a user drops their own
    complete .mol2 structures in there (or points this at any other
    directory, e.g. a personal collection kept outside the repo).

    Returns a list of (category_name, [(display_name, file_path), ...])
    tuples, both levels sorted alphabetically, skipping any category folder
    that turns out to contain no .mol2 files at all. """
    if base_dir is None:
        base_dir = os.path.join ( os.path.dirname ( os.path.abspath ( __file__ ) ), "structures" )

    categories = [ ]
    if not os.path.isdir ( base_dir ):
        return categories

    for category_name in sorted ( os.listdir ( base_dir ) ):
        category_path = os.path.join ( base_dir, category_name )
        if not os.path.isdir ( category_path ):
            continue
        entries = [ ]
        for filename in sorted ( os.listdir ( category_path ) ):
            if filename.lower ( ).endswith ( ".mol2" ):
                display_name = os.path.splitext ( filename )[0]
                entries.append ( ( display_name, os.path.join ( category_path, filename ) ) )
        if entries:
            categories.append ( ( category_name, entries ) )

    return categories


def _resolve_symbol ( atom_name, atom_type ):
    """ [EN] Deliberately OPPOSITE column priority from fragment_library.
    _resolve_symbol() (which tries the "type" column FIRST): here the
    "name" column (trailing digits stripped, e.g. "Cl1" -> "Cl") is tried
    FIRST, falling back to the type column's pre-dot part (e.g. "O.3" ->
    "O") only if the name doesn't resolve.

    [EN] BUG FOUND live-testing against this checkout's own examples/mol2/
    valacyclovir.mol2 (a real, complete-molecule MOL2 export from Gabedit,
    AMBER/GAFF-typed): type-column-first, matched against this module's
    OWN broader _PERIODIC_TABLE_SYMBOLS (every element, unlike fragment_
    library.py's small hand-picked subset), silently mis-resolved several
    atoms -- AMBER atom TYPE codes "NA"/"CA"/"OS"/"NB" (nitrogen-aromatic/
    alpha-carbon-style/ester-oxygen/aromatic-nitrogen-b, in AMBER's own
    naming, nothing to do with real elements) happen to be VALID element
    symbols too (Sodium/Calcium/Osmium/Niobium) -- fragment_library.py's
    own small subset never had this exposure (none of those 4 elements are
    in its supported set at all), but this module's intentionally broader
    table does collide. This file's plain element-letter NAME column
    ("N"/"C"/"O"/"H", no residue-role naming) resolves every atom
    correctly, so trying it first fixes this exact file.

    [EN] KNOWN, DOCUMENTED LIMITATION (not fixed this pass, no concrete
    failing file to verify a fix against): a structure file using
    PROTEIN-BACKBONE-STYLE atom NAMES (e.g. "CA" meaning alpha-carbon,
    "CB" beta-carbon -- common for amino-acid/residue templates, NOT for
    the small-molecule/ligand/cofactor files this feature's own "biblioteca
    de estruturas prontas" is primarily meant for) would hit the exact same
    class of collision from the OTHER column instead (name "CA" -> wrongly
    resolves to Calcium) -- unavoidable with this pure structural/textual
    heuristic (no dictionary of standard residue atom names is consulted).
    If that ever becomes a real, hit case, the fix is a residue-name-aware
    exception list, not another column-priority flip. """
    for candidate in ( atom_name.rstrip ( "0123456789" ), atom_type.split ( "." )[0] ):
        resolved = _PERIODIC_TABLE_SYMBOLS.get ( candidate.strip ( ).upper ( ) )
        if resolved is not None:
            return resolved
    return None


def _parse_mol2 ( path ):
    """ [EN] Same minimal Tripos MOL2 reader as fragment_library._parse_
    mol2() (only the 3 sections a structure file needs: MOLECULE/ATOM/
    BOND -- no charges, no crystal symmetry, no sequence/residue handling),
    kept as its own copy rather than imported from fragment_library.py so
    this module has zero dependency on that one (element resolution differs
    -- see _resolve_symbol() above -- and the two are conceptually
    independent features that just happen to share a file format).

    Returns ( atoms, bonds ):
      atoms -- list of (symbol, x, y, z), 0-indexed, in file order.
      bonds -- list of (i, j, mol2_bond_code), 0-indexed atom positions,
               mol2_bond_code the RAW string from the file (e.g. "1", "ar").

    Raises StructureError for any structural problem (missing sections,
    wrong atom/bond counts, unresolvable element). """
    with open ( path, "r" ) as handle:
        lines = [ line.rstrip ( "\n" ) for line in handle ]

    n_atoms = None
    n_bonds = None
    atoms   = [ ]
    bonds   = [ ]

    i = 0
    while i < len ( lines ):
        line = lines[i].strip ( )

        if line == "@<TRIPOS>MOLECULE":
            counts = lines[i + 2].split ( )
            n_atoms = int ( counts[0] )
            n_bonds = int ( counts[1] ) if len ( counts ) > 1 else 0
            i += 3
            continue

        elif line == "@<TRIPOS>ATOM":
            if n_atoms is None:
                raise StructureError ( "{}: ATOM section appears before MOLECULE section.".format ( path ) )
            for k in range ( n_atoms ):
                i += 1
                fields = lines[i].split ( )
                if len ( fields ) < 6:
                    raise StructureError ( "{}: malformed ATOM line {!r}.".format ( path, lines[i] ) )
                atom_name = fields[1]
                x, y, z   = float ( fields[2] ), float ( fields[3] ), float ( fields[4] )
                atom_type = fields[5]
                symbol = _resolve_symbol ( atom_name, atom_type )
                if symbol is None:
                    raise StructureError ( "{}: atom #{} ({!r}/{!r}) is not a recognised element.".format (
                            path, k + 1, atom_name, atom_type ) )
                atoms.append ( ( symbol, x, y, z ) )
            i += 1
            continue

        elif line == "@<TRIPOS>BOND":
            for k in range ( n_bonds ):
                i += 1
                fields = lines[i].split ( )
                if len ( fields ) < 4:
                    raise StructureError ( "{}: malformed BOND line {!r}.".format ( path, lines[i] ) )
                atom_i = int ( fields[1] ) - 1
                atom_j = int ( fields[2] ) - 1
                bonds.append ( ( atom_i, atom_j, fields[3] ) )
            i += 1
            continue

        else:
            i += 1

    if n_atoms is None:
        raise StructureError ( "{}: no @<TRIPOS>MOLECULE section found.".format ( path ) )
    if len ( atoms ) != n_atoms:
        raise StructureError ( "{}: expected {} atoms, found {}.".format ( path, n_atoms, len ( atoms ) ) )
    if len ( bonds ) != n_bonds:
        raise StructureError ( "{}: expected {} bonds, found {}.".format ( path, n_bonds, len ( bonds ) ) )

    return atoms, bonds


def load_structure ( path ):
    """ [EN] Parses a COMPLETE-molecule .mol2 file, returning a plain dict
    ready for atom_ops.add_structure_at_position():

        { "name"  : str  (the file's own basename, no extension),
          "atoms" : [ (symbol, x, y, z), ... ]              0-indexed,
          "bonds" : [ (i, j, order, is_aromatic), ... ]     0-indexed,
                    order already resolved to a plain int (1/2/3) }

    Unlike fragment_library.load_fragment(), there is NO attachment-point
    requirement/validation at all -- this loader trusts the file's own
    atoms/bonds as a complete, ready-to-place structure. Only structural
    problems (malformed sections, an unresolvable element, an unsupported
    bond code) raise StructureError; a fully-capped molecule with every
    valence already satisfied -- the normal, expected case here -- is
    exactly what this loads without complaint. """
    atoms, raw_bonds = _parse_mol2 ( path )

    if not atoms:
        raise StructureError ( "{}: structure has no atoms.".format ( path ) )

    bonds = [ ]
    for ( i, j, code ) in raw_bonds:
        normalized_code = code.strip ( ).lower ( )
        order = _MOL2_BOND_ORDERS.get ( normalized_code )
        if order is None:
            raise StructureError ( "{}: bond {}-{} has unsupported type {!r}.".format (
                    path, i + 1, j + 1, code ) )
        bonds.append ( ( i, j, order, normalized_code in _MOL2_AROMATIC_CODES ) )

    return {
        "name"  : os.path.splitext ( os.path.basename ( path ) )[0],
        "atoms" : atoms,
        "bonds" : bonds,
    }
