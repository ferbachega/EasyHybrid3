#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  fragment_library.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Builder fragment library -- small, reusable pieces (functional groups for
#  now; rings/residues are later phases, see the plan file) that can be
#  attached onto a hydrogen atom of the molecule currently being edited
#  (see atom_ops.attach_fragment_at_hydrogen() / click_mode.handle_click_to_
#  attach_fragment()).
#
#  [EN] FRAGMENT FILE CONVENTION (the whole point of this module -- what
#  makes the library "easy to customize", per the user's own request): a
#  valid attachable fragment .mol2 file is just a normal, chemically
#  COMPLETE small molecule EXCEPT for exactly ONE atom left one bond short
#  of its standard valence (e.g. -OH is written as a plain O-H, -NH2 as
#  N(H)(H), -COOH as the complete carboxylic acid with the carbonyl
#  carbon's 4th bond left open). That single under-valent atom is
#  auto-detected here as the attachment point -- no dummy/placeholder
#  atoms, no filename or atom-naming convention to remember: take any
#  complete molecule, delete one terminal hydrogen, save as .mol2.
#
#  [EN] Deliberately does NOT reuse util/extras/MOL2FileReader.py: that
#  reader is unreachable dead code in this checkout -- src/util/extras/
#  __init__.py imports ChargeRestraintModel/CPHFSolver/DFTDefinitions/etc,
#  files that belong to an unrelated pDynamo QC-model package and simply
#  don't exist in src/util/extras/ (confirmed via `ls`), so `import
#  util.extras` (needed even for a direct submodule import, since
#  MOL2FileReader.py itself does a relative `from .ExportImport import
#  _Importer`) raises ModuleNotFoundError every time -- confirmed by
#  actually running it, not just reading the source. This explains why
#  grep finds zero real call sites for that reader anywhere in the live
#  app despite its own _Importer.AddHandler(...) registration: nothing can
#  actually reach it. Writing a small, self-contained parser here instead
#  avoids that landmine entirely and only needs to understand the 3 MOL2
#  sections a fragment file actually uses (MOLECULE/ATOM/BOND) -- no
#  charges, no crystal symmetry, no sequence/residue handling.
# ============================================================================
import os


# [EN] Maps the UPPER-CASE form of every element this Builder's own
# atom_ops.STANDARD_VALENCE table supports to its PROPER-CASE symbol (e.g.
# "CL" -> "Cl") -- used both to resolve a MOL2 atom's element and to decide
# which atoms can even be valence-checked at all (anything outside this set
# is rejected with a clear error rather than silently guessed at, matching
# atom_ops.adjust_hydrogens()'s own "unknown elements are left alone rather
# than guessed at" philosophy -- except here, for a FRAGMENT DEFINITION,
# "left alone" isn't safe, so it's a hard error instead). Kept as its own
# copy (not imported from atom_ops.STANDARD_VALENCE directly) so this
# module has zero import-time dependency on atom_ops -- only load_fragment()
# needs it, and importing atom_ops just for a dict would pull in this
# entire Builder's GTK/pDynamo dependency chain for no real benefit.
_SUPPORTED_ELEMENTS = { 'H': 'H', 'C': 'C', 'N': 'N', 'O': 'O', 'F': 'F',
                         'CL': 'Cl', 'BR': 'Br', 'I': 'I', 'P': 'P', 'S': 'S' }

_STANDARD_VALENCE = { 'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1,
                       'CL': 1, 'BR': 1, 'I': 1, 'P': 3, 'S': 2 }

# [EN] Tripos MOL2 bond-type codes this Builder can represent -- "am"
# (amide) is chemically just a single bond (the code only exists so
# force-field tools can flag amide conjugation separately; this Builder has
# no such concept, so it's treated as a plain single bond). "ar"
# (aromatic), "du"/"un" (dummy/unknown) and anything else are explicitly
# NOT supported in this pass -- see this module's own top-of-file note on
# why aromatic/ring fragments are out of scope (no Kekulization step exists
# to turn "ar" into a concrete Single/Double order).
_MOL2_BOND_ORDERS = { '1': 1, '2': 2, '3': 3, 'am': 1 }


class FragmentError ( Exception ):
    """ Raised by load_fragment() for any malformed/unsupported fragment
    file -- always carries a short, specific, user-facing message (which
    atom, which problem) rather than a generic parse failure, since the
    whole point of this module is helping someone author their OWN
    fragment files correctly. """
    pass


def discover_fragment_files ( base_dir = None ):
    """ [EN] Walks base_dir's IMMEDIATE subdirectories (each one a
    category -- "functional_groups" for this pass, "rings"/"residues" in
    later phases, per the plan file) and returns every ".mol2" file found
    directly inside each one (not recursing further -- a category is a
    flat list of fragments, not its own tree).

    base_dir defaults to this module's own "fragments" subdirectory, so
    the shipped starter set is found with no configuration -- but a
    caller can point this at any other directory too (e.g. a user's own
    personal fragment collection kept outside the repo), which is the
    other half of "easy to customize" alongside the file-convention itself.

    Returns a list of (category_name, [(display_name, file_path), ...])
    tuples, both levels sorted alphabetically, skipping any category
    folder that turns out to contain no .mol2 files at all. """
    if base_dir is None:
        base_dir = os.path.join ( os.path.dirname ( os.path.abspath ( __file__ ) ), "fragments" )

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
    """ [EN] Resolves a MOL2 atom's element symbol from its "type" column
    (SYBYL convention, e.g. "O.3"/"C.ar"/"N.1" -- take the part before the
    dot) first, falling back to the "name" column with any trailing
    digits stripped (e.g. "Cl1" -> "Cl", "H12" -> "H") if the type column
    didn't resolve to a known element. Returns the PROPER-CASE symbol
    (e.g. "Cl", not "CL"/"cl"), or None if neither column resolves to an
    element this Builder supports. """
    for candidate in ( atom_type.split ( "." )[0], atom_name.rstrip ( "0123456789" ) ):
        resolved = _SUPPORTED_ELEMENTS.get ( candidate.strip ( ).upper ( ) )
        if resolved is not None:
            return resolved
    return None


def _parse_mol2 ( path ):
    """ [EN] Minimal Tripos MOL2 reader -- only the 3 sections a fragment
    file actually needs (see this module's own top-of-file note for why
    this exists instead of reusing util/extras/MOL2FileReader.py).

    Returns ( atoms, bonds ):
      atoms -- list of (symbol, x, y, z), 0-indexed, in file order.
      bonds -- list of (i, j, mol2_bond_code), 0-indexed atom positions,
               mol2_bond_code the RAW string from the file (e.g. "1",
               "ar") -- order/validity is resolved later by the caller,
               not here, so a malformed bond code can be reported with
               the fragment's own name/context instead of a bare parse
               error.

    Raises FragmentError for any structural problem (missing sections,
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
            # [EN] line i+1 = molecule name (unused here, just skipped);
            # line i+2 = "n_atoms n_bonds" (n_bonds may be omitted for a
            # single-atom, zero-bond fragment -- see e.g. this library's
            # own chloro.mol2/bromo.mol2/fluoro.mol2 starter fragments).
            counts = lines[i + 2].split ( )
            n_atoms = int ( counts[0] )
            n_bonds = int ( counts[1] ) if len ( counts ) > 1 else 0
            i += 3
            continue

        elif line == "@<TRIPOS>ATOM":
            if n_atoms is None:
                raise FragmentError ( "{}: ATOM section appears before MOLECULE section.".format ( path ) )
            for k in range ( n_atoms ):
                i += 1
                fields = lines[i].split ( )
                if len ( fields ) < 6:
                    raise FragmentError ( "{}: malformed ATOM line {!r}.".format ( path, lines[i] ) )
                atom_name = fields[1]
                x, y, z   = float ( fields[2] ), float ( fields[3] ), float ( fields[4] )
                atom_type = fields[5]
                symbol = _resolve_symbol ( atom_name, atom_type )
                if symbol is None:
                    raise FragmentError ( "{}: atom #{} ({!r}/{!r}) is not a supported element.".format (
                            path, k + 1, atom_name, atom_type ) )
                atoms.append ( ( symbol, x, y, z ) )
            i += 1
            continue

        elif line == "@<TRIPOS>BOND":
            for k in range ( n_bonds ):
                i += 1
                fields = lines[i].split ( )
                if len ( fields ) < 4:
                    raise FragmentError ( "{}: malformed BOND line {!r}.".format ( path, lines[i] ) )
                atom_i = int ( fields[1] ) - 1
                atom_j = int ( fields[2] ) - 1
                bonds.append ( ( atom_i, atom_j, fields[3] ) )
            i += 1
            continue

        else:
            i += 1

    if n_atoms is None:
        raise FragmentError ( "{}: no @<TRIPOS>MOLECULE section found.".format ( path ) )
    if len ( atoms ) != n_atoms:
        raise FragmentError ( "{}: expected {} atoms, found {}.".format ( path, n_atoms, len ( atoms ) ) )
    if len ( bonds ) != n_bonds:
        raise FragmentError ( "{}: expected {} bonds, found {}.".format ( path, n_bonds, len ( bonds ) ) )

    return atoms, bonds


def load_fragment ( path ):
    """ [EN] Parses and VALIDATES a fragment .mol2 file, returning a
    plain dict ready for atom_ops.attach_fragment_at_hydrogen():

        { "name"      : str  (the file's own basename, no extension),
          "atoms"     : [ (symbol, x, y, z), ... ]   0-indexed,
          "bonds"     : [ (i, j, order), ... ]        0-indexed, order is
                        already resolved to a plain int (1/2/3),
          "root_index": int   the auto-detected attachment-point atom }

    Validation (each failure raises FragmentError with a specific,
    actionable message -- see this module's own top-of-file note on why
    that matters for a "make your own fragments" feature):
      - every atom must resolve to an element in atom_ops.STANDARD_VALENCE
        (already enforced by _parse_mol2()/_resolve_symbol()).
      - every bond code must be a supported, non-aromatic type (_MOL2_
        BOND_ORDERS) -- aromatic/dummy/unknown bonds are rejected outright,
        not approximated.
      - EXACTLY ONE atom may be under-valent, and by EXACTLY one bond
        order unit (this Builder only ever attaches via a single bond,
        see attach_fragment_at_hydrogen()'s own docstring) -- zero
        under-valent atoms (a fully-capped molecule, no attachment point)
        or more than one (multi-point attachment, not supported this
        pass) are both rejected. """
    atoms, raw_bonds = _parse_mol2 ( path )

    bonds = [ ]
    for ( i, j, code ) in raw_bonds:
        order = _MOL2_BOND_ORDERS.get ( code )
        if order is None:
            raise FragmentError ( "{}: bond {}-{} has unsupported type {!r} (aromatic/ring fragments "
                                   "aren't supported yet -- see the roadmap).".format (
                                   path, i + 1, j + 1, code ) )
        bonds.append ( ( i, j, order ) )

    bond_order_sum = [ 0 ] * len ( atoms )
    for ( i, j, order ) in bonds:
        bond_order_sum[i] += order
        bond_order_sum[j] += order

    open_valence_atoms = [ ]
    for idx, ( symbol, _x, _y, _z ) in enumerate ( atoms ):
        target = _STANDARD_VALENCE[ symbol.upper ( ) ]
        shortfall = target - bond_order_sum[idx]
        if shortfall > 0:
            open_valence_atoms.append ( ( idx, shortfall ) )

    if not open_valence_atoms:
        raise FragmentError ( "{}: fragment has no open attachment point (every atom's valence is "
                               "already satisfied) -- remove one hydrogen from the intended "
                               "attachment atom.".format ( path ) )
    if len ( open_valence_atoms ) > 1:
        raise FragmentError ( "{}: fragment has {} open attachment points (expected exactly 1) -- "
                               "multi-point attachment isn't supported yet.".format (
                               path, len ( open_valence_atoms ) ) )

    root_index, shortfall = open_valence_atoms[0]
    if shortfall != 1:
        raise FragmentError ( "{}: the attachment atom is {} bonds short of its standard valence "
                               "(expected exactly 1) -- this Builder only attaches fragments via a "
                               "single bond.".format ( path, shortfall ) )

    return {
        "name"      : os.path.splitext ( os.path.basename ( path ) )[0],
        "atoms"     : atoms,
        "bonds"     : bonds,
        "root_index": root_index,
    }
