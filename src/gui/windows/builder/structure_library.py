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
import re
from vismol.utils.elements import PeriodicTable
from vismol.utils.debug import dprint


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


# . 2026-10-24: readable names for the library trees ("01_aromatic_carbocyclic_rings" ->
#   "Aromatic carbocyclic rings", "barbituric_acid" -> "barbituric acid"); the folders and
#   files keep their names.
_NAME_FIXES = { "sulfurvi": "sulfur(VI)", "3d": "3D", "4d": "4d", "5d": "5d", "n": "N", "o": "O", "s": "S" }

def pretty_library_name ( name, category = False ):
    import re as _re
    text = _re.sub ( r"^\d+_", "", str ( name ) )
    text = _re.sub ( r"(^|_)o_s(_|$)", r"\1O/S\2", text ).replace ( "6membered", "6-membered" )
    words = [ ]
    for w in text.split ( "_" ):
        if not w:
            continue
        if w in ( "3d", "4d", "5d" ) and words and words[-1] == "metals":
            words.append ( w )
        else:
            words.append ( _NAME_FIXES.get ( w, w ) )
    text = " ".join ( words )
    if category and text:
        text = text[0].upper ( ) + text[1:]
    return text


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
    exception list, not another column-priority flip.

    [EN] 2026-10-17 BUG FIX (user's report: ethanolamine.mol2's "CA" atom,
    type "C.3", came in as Calcium -- so DYFF + xTB ran on a calcium
    complex). The SYBYL type column IS unambiguous in two cases, so it
    now wins there, before the name: a dotted type ("C.3", "N.ar", "O.co2"
    -- AMBER/GAFF types never have a dot) and a type written exactly as
    an element symbol ("C", "H", "Cl", "Ca", "Zn" -- AMBER types are
    upper case, "CA"/"NA"/"OS", GAFF lower case, "ca"/"c3"). Only the
    other types (AMBER/GAFF) fall back to name-first, as before. """
    base_type = atom_type.strip ( ).split ( "." )[0]
    resolved = _PERIODIC_TABLE_SYMBOLS.get ( base_type.upper ( ) )
    if resolved is not None and ( "." in atom_type or base_type == resolved ):
        return resolved
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


def _parse_formal_charges ( path ):
    """ [EN] 2026-10-13: { 0-based atom index: formal charge } from the file's
    @<TRIPOS>UNITY_ATOM_ATTR section (the standard Tripos place for FORMAL
    charges, same convention as fragment_library: blocks of "<atom_id>
    <n_attrs>" followed by "<name> <value>" lines, only "charge" used). Added
    for the transition-metal complexes (category 31_*), whose donor and metal
    formal charges make DYFF give the metal its oxidation state (covalent
    model: X donors neutral, L donors +1, the metal the rest). A file
    without the section -> { } (every atom neutral, as before). """
    charges = { }
    with open ( path, "r" ) as handle:
        lines = [ line.rstrip ( "\n" ) for line in handle ]
    i = 0
    while i < len ( lines ):
        if lines[i].strip ( ) == "@<TRIPOS>UNITY_ATOM_ATTR":
            i += 1
            while i < len ( lines ) and not lines[i].strip ( ).startswith ( "@<TRIPOS>" ):
                fields = lines[i].split ( )
                if len ( fields ) == 2 and fields[0].isdigit ( ) and fields[1].isdigit ( ):
                    atom_index = int ( fields[0] ) - 1
                    for _ in range ( int ( fields[1] ) ):
                        i += 1
                        attr = lines[i].split ( ) if i < len ( lines ) else [ ]
                        if len ( attr ) == 2 and attr[0].lower ( ) == "charge":
                            try:
                                charges[atom_index] = int ( float ( attr[1] ) )
                            except ValueError:
                                raise StructureError ( "{}: bad formal charge {!r} for atom #{}.".format ( path, attr[1], atom_index + 1 ) )
                i += 1
            continue
        i += 1
    return charges


# . target valence (sum of bond orders) of an atom with formal charge q, for the Kekulization
_KEKULE_VALENCE = { "C": lambda q: 4 - abs ( q ), "N": lambda q: 3 + q, "O": lambda q: 2 + q, "S": lambda q: 2 + q,
                    "Se": lambda q: 2 + q, "P": lambda q: 3 + q, "B": lambda q: 3 + q }

def _kekulize_aromatic ( symbols, bonds, formal ):
    """ Orders for every bond (list, same order as bonds [ i, j, order, aromatic ])
        with the aromatic bonds turned into single/double so that every atom that
        lacks valence gets exactly one double bond -- explicit orders and formal
        charges respected (maximum matching on the aromatic bonds). None when no
        such assignment exists (or networkx is missing). """
    try:
        import networkx as nx
    except ImportError:
        return None
    used = [ 0 ] * len ( symbols )
    for i, j, order, aromatic in bonds:
        used[i] += order; used[j] += order
    deficit = [ ]
    for a, symbol in enumerate ( symbols ):
        rule = _KEKULE_VALENCE.get ( symbol )
        deficit.append ( rule ( formal[a] ) - used[a] if rule is not None else 0 )
    graph = nx.Graph ( )
    needing = set ( )
    for k, ( i, j, order, aromatic ) in enumerate ( bonds ):
        if aromatic and deficit[i] >= 1 and deficit[j] >= 1:
            graph.add_edge ( i, j, bond = k ); needing.update ( ( i, j ) )
    for i, j, order, aromatic in bonds:
        for a in ( i, j ):
            if aromatic and deficit[a] >= 1: needing.add ( a )
    if any ( deficit[a] > 1 for a in needing ):
        return None
    matching = nx.max_weight_matching ( graph, maxcardinality = True )
    matched = { a for pair in matching for a in pair }
    if matched != needing:
        return None
    orders = [ order for ( i, j, order, aromatic ) in bonds ]
    for a, b in matching:
        orders[graph[a][b]["bond"]] = 2
    return orders


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
        bonds.append ( [ i, j, order, normalized_code in _MOL2_AROMATIC_CODES ] )

    # [EN] REAL BUG FOUND AND FIXED, 2026-09-24 -- user's own report: "as
    # ligacoes nao estao devidamente atribuidas, todas estao como simples,
    # desrespeitando a ordem de ligacao do mol2". Traced end to end
    # (parsing -> atom_ops.add_structure_at_position() -> manual_bond_
    # orders -> rendering): explicit "1"/"2"/"3" mol2 bond codes were
    # ALREADY passing through correctly (verified against real files:
    # nitrile.mol2's C#N "3", urea.mol2's C=O "2"). The actual cause: "ar"
    # (aromatic) bonds are, BY DESIGN (see this module's own top-of-file
    # note), stored as order=1 + a separate aromatic flag -- mirroring the
    # Builder's own manual "Aromatic" sidebar tool -- but nothing ever
    # Kekulizes that flag into a concrete alternating single/double
    # pattern for RENDERING (SticksRepresentation only ever reads `bond.
    # bond_order`, never the aromatic flag -- see representations.py).
    # Since most of this library's own shipped content (benzene rings,
    # pyridines, purines, fused medchem scaffolds...) is essentially 100%
    # "ar"-coded on their ring bonds, placing one of THOSE structures
    # really did show "every bond single", exactly matching the report --
    # while an explicit-double/triple structure (urea, nitrile) already
    # rendered correctly, which is why this looked inconsistent rather
    # than uniformly broken.
    #
    # Fixed by reusing this app's OWN already-proven Wang&Case bond-order
    # perceiver (vismol.core.bond_order_perception.perceive_bond_orders --
    # the exact same one atom_ops.bootstrap_manual_bonds_from_existing()
    # already uses to Kekulize a whole "Edit in Builder" molecule from
    # scratch) over this structure's OWN local atoms/bonds, but ONLY to
    # resolve the "ar"-coded ones -- explicit "1"/"2"/"3" bonds from the
    # file are left untouched (already correct, no reason to second-guess
    # the file's own explicit data), while still feeding the perceiver the
    # WHOLE local bond graph (not just the aromatic subset) so its own
    # valence-penalty minimisation has full context to alternate correctly
    # around the ring. Same "never let a perception failure block loading"
    # guard as bootstrap_manual_bonds_from_existing() -- falls back to the
    # OLD order=1 behaviour (not a crash) if perception itself raises.
    charges = _parse_formal_charges ( path )
    for k in charges:
        if not 0 <= k < len ( atoms ):
            raise StructureError ( "{}: formal charge given for non-existent atom #{}.".format ( path, k + 1 ) )
    formal = [ charges.get ( k, 0 ) for k in range ( len ( atoms ) ) ]

    # [EN] 2026-10-21: the "ar" bonds are first Kekulized with the file's explicit
    # orders and formal charges as constraints (_kekulize_aromatic). The perceiver
    # below ignores both: on files that mix explicit "2" and "ar" bonds (porphine,
    # protoporphyrin, cobalamins, F420) it put a second double bond on carbons that
    # already had one (C with 5 bonds). The perceiver stays as the fallback.
    kekule = _kekulize_aromatic ( [ atom[0] for atom in atoms ], bonds, formal ) \
             if any ( is_aromatic for ( _i, _j, _order, is_aromatic ) in bonds ) else None
    if kekule is not None:
        for entry, order in zip ( bonds, kekule ):
            entry[2] = order
    elif any ( is_aromatic for ( _i, _j, _order, is_aromatic ) in bonds ):
        try:
            from vismol.core.bond_order_perception import perceive_bond_orders
            symbols = [ atom[0] for atom in atoms ]
            pairs = [ ( i, j ) for ( i, j, _order, _arom ) in bonds ]
            order_map, _tps = perceive_bond_orders ( symbols, pairs )
            for entry in bonds:
                i, j, order, is_aromatic = entry
                if is_aromatic:
                    key = ( i, j ) if i < j else ( j, i )
                    entry[2] = int ( order_map.get ( key, order ) )
        except Exception as exc:
            dprint ( "WARNING structure_library.load_structure: bond order perception "
                    "failed for {} ({}) -- aromatic bonds kept as plain single order.".format ( path, exc ) )

    bonds = [ tuple ( entry ) for entry in bonds ]

    return {
        "name"  : os.path.splitext ( os.path.basename ( path ) )[0],
        "atoms" : atoms,
        "bonds" : bonds,
        "formal_charges" : formal,
    }


#=============================================================================
#  Charge and spin of a library structure (2026-10-21)
#=============================================================================
#  [EN] User's request: the Structures tab shows, for the selected structure, its
#  total charge and spin multiplicity -- above all for the transition-metal
#  complexes. The charge is the sum of the file's formal charges. The
#  multiplicity: structures/structure_info.json when the structure has an
#  entry (metal cofactors, curated); the transition-metal complexes of the
#  31-33 folders come from tmQM, whose structures are all closed-shell DFT
#  references (multiplicity 1) -- their metal and oxidation state are in the
#  file name ("Fe_II_octahedral_PIMMOQ"); any other molecule: 1 for an even
#  number of electrons, 2 for an odd one.

_ROMAN = { "0": 0, "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7, "VIII": 8 }
_D_BLOCK = { s: g for row in ( ( "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn" ),
                               ( "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd" ),
                               ( "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg" ) )
             for g, s in zip ( range ( 3, 13 ), row ) }
_INFO_CACHE = { }

def _structure_info_table ( base_dir ):
    path = os.path.join ( base_dir, "structure_info.json" )
    if path not in _INFO_CACHE:
        import json
        try:
            with open ( path ) as handle:
                _INFO_CACHE[path] = { k: v for k, v in json.load ( handle ).items ( ) if not k.startswith ( "_" ) }
        except ( OSError, ValueError ):
            _INFO_CACHE[path] = { }
    return _INFO_CACHE[path]

def _d_electrons ( metal, oxidation ):
    group = _D_BLOCK.get ( metal )
    return None if group is None or oxidation is None else group - oxidation

def describe_structure ( structure, path ):
    """ { "charge", "multiplicity", "source", "metal", "oxidation", "d_electrons",
          "high_spin", "note" } for a structure returned by load_structure(). """
    from vismol.utils.elements import PeriodicTable as _PT
    table = _PT ( ).elements_by_symbol
    charge = int ( sum ( structure.get ( "formal_charges" ) or [ ] ) )
    electrons = sum ( int ( table[a[0]][0] ) for a in structure["atoms"] if a[0] in table ) - charge
    folder = os.path.basename ( os.path.dirname ( path ) )
    key = "{}/{}".format ( folder, os.path.basename ( path ) )
    entry = _structure_info_table ( os.path.dirname ( os.path.dirname ( path ) ) ).get ( key, { } )
    metals = sorted ( { a[0] for a in structure["atoms"] if a[0] in _D_BLOCK } )
    info = { "charge": charge, "multiplicity": 2 if electrons % 2 else 1, "source": "electron count",
             "metal": metals[0] if len ( metals ) == 1 else ( ", ".join ( metals ) or None ),
             "oxidation": None, "d_electrons": None, "high_spin": None, "note": None }
    match = re.match ( r"^([A-Z][a-z]?)_(0|[IVX]+)_", os.path.basename ( path ) )
    if match and match.group ( 1 ) in _D_BLOCK and folder[:2] in ( "31", "32", "33" ):
        info.update ( metal = match.group ( 1 ), oxidation = _ROMAN.get ( match.group ( 2 ) ), multiplicity = 1,
                      source = "tmQM reference (closed-shell DFT)" )
    if entry:
        info.update ( { k: entry[k] for k in ( "metal", "oxidation", "multiplicity", "note" ) if k in entry } )
        info["source"] = "curated"
    if info["metal"] in _D_BLOCK and info["oxidation"] is not None:
        d = _d_electrons ( info["metal"], info["oxidation"] )
        info["d_electrons"] = d
        # . high spin only matters for 3d metals; 4d/5d and organometallic eta
        #   (18-electron) complexes are low spin
        three_d = _D_BLOCK.get ( info["metal"] ) is not None and info["metal"] in ( "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn" )
        if d is not None and 2 <= d <= 8 and three_d and "eta_complex" not in os.path.basename ( path ):
            info["high_spin"] = ( d if d <= 5 else 10 - d ) + 1
    return info

_MULTIPLICITY_NAMES = { 1: "singlet", 2: "doublet", 3: "triplet", 4: "quartet", 5: "quintet", 6: "sextet" }

def describe_structure_text ( structure, path ):
    """ The text of the Structures tab's status label. """
    info = describe_structure ( structure, path )
    m = info["multiplicity"]
    lines = [ "{} · {} atoms".format ( structure["name"], len ( structure["atoms"] ) ),
              "Charge {:+d} · multiplicity {} ({})".format ( info["charge"], m, _MULTIPLICITY_NAMES.get ( m, "" ) ).replace ( "+0", "0" ) ]
    if info["metal"] and info["oxidation"] is not None:
        roman = { v: k for k, v in _ROMAN.items ( ) }.get ( info["oxidation"], str ( info["oxidation"] ) )
        text = "{}({})".format ( info["metal"], roman )
        if info["d_electrons"] is not None: text += ", d{}".format ( info["d_electrons"] )
        if info["source"].startswith ( "tmQM" ):
            text += " · multiplicity from the tmQM reference (closed shell)"
        if info["high_spin"] and info["high_spin"] != m:
            text += " · high spin would be {}".format ( info["high_spin"] )
        lines.append ( text )
    if info["note"]:
        lines.append ( info["note"] )
    return "\n".join ( lines )
