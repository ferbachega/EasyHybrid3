#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: "Add Missing Hydrogens" -- residue-template-driven H detection/building
#
""" Finds and adds missing hydrogen atoms to standard amino-acid residues
    of an already-loaded pDynamo System, reusing pDynamo3's OWN chemistry
    data/machinery instead of re-implementing functional-group detection
    or geometry building from scratch:

    - pBabel.PDBComponentLibrary supplies, per residue name, the COMPLETE
      expected atom list (every heavy atom AND every named hydrogen) plus
      named protonation "variants" (e.g. HIS's "Delta Protonated" /
      "Epsilon Protonated" / "Doubly Protonated" -- HID/HIE/HIP; ASP/GLU's
      "Deprotonated" default). Missing atoms are found by DIFFING this
      real template against a residue's actual atom names -- no custom
      carbon/carbonyl/hydroxyl/amide pattern matching needed.
    - pSimulation.BuildHydrogenCoordinates3FromConnectivity builds 3D
      coordinates for any atom flagged "undefined" using only the bond
      graph + per-element coordination geometry -- no forcefield/mmModel
      required (confirmed: this module's own docstring states energy
      models are not needed).

    pDynamo3 ships no CYM (deprotonated Cys) / LYN (neutral Lys) variant
    -- CUSTOM_VARIANTS below adds both, in the exact same shape as the
    library's own shipped ASP/GLU "Deprotonated" variants
    (pBabel/PDBComponentLibrary.py:206-207): one heavy-atom-bonded H
    removed, one formal charge adjusted.

    SCOPE (deliberate, v1): only HYDROGEN atoms are ever added. "H2",
    "OXT" and "HXT" are always excluded from a residue's expected atom
    set -- these three names only apply to a real chain terminus (via
    pDynamo3's own "Peptide"/"N Terminal"/"C Terminal" link/variant
    machinery, layered on top of whatever protonation variant is chosen),
    and OXT is a heavy atom besides. Capping real N-/C-termini is a
    separate, narrower concern this tool does not attempt -- for every
    ordinary (interior) residue, excluding these three names exactly
    matches what pDynamo3's own "Peptide" link would do anyway, so this
    is not a shortcut on correctness for the overwhelming majority of
    real residues, just an explicit, reported limitation for whichever
    residue happens to be the very first/last of its chain.

    Import order note (this project's own known pDynamo3 quirk): pBabel
    must be imported before pMolecule or a circular ImportError results
    ("cannot import name 'EnergyModelPriority'"). All pDynamo imports
    here are lazy (inside functions), so this module stays importable
    standalone -- matching vina_runner.py / vismol_pdbqt_builder.py's own
    convention -- and relies on the caller (EasyHybrid itself) having
    already imported pDynamo3 in the right order by the time these
    functions actually run.

    NON-STANDARD residues (ligands, cofactors, any organic molecule with
    no PDBComponentLibrary entry) are handled by a SEPARATE path,
    protonate_non_standard_residue(), that shells out to OpenBabel (the
    same external tool -- and same `-p <ph>` flag -- already used for
    ligand docking prep in vina_runner.py) instead of a residue-name
    template: OpenBabel perceives valence AND pH-dependent protonation
    state for an arbitrary organic molecule from geometry alone, which
    a name-keyed template fundamentally cannot do for something with no
    name in the library.
"""

import os
import shutil
import subprocess
import tempfile

# . Atom names that only ever belong to a genuine chain terminus -- see
# the module docstring's SCOPE note. Never added by this tool.
_TERMINUS_ONLY_ATOM_NAMES = frozenset(('H2', 'OXT', 'HXT'))

# . Residues with a real protonation choice worth surfacing in a GUI.
# 'Protonated' (or any label not naming an actual variant) means "use
# the component's own unmodified/base atom list" -- represented
# internally as variant_label=None.
PROTONATION_VARIANTS = {
    'HIS': ['Delta Protonated', 'Epsilon Protonated', 'Doubly Protonated'],
    'ASP': ['Protonated', 'Deprotonated'],
    'GLU': ['Protonated', 'Deprotonated'],
    'CYS': ['Protonated', 'Deprotonated'],
    'LYS': ['Protonated', 'Neutral'],
}

_BASE_VARIANT_LABELS = frozenset(('Protonated',))

_LIBRARY = None


def _get_library():
    """ Lazily builds (once per process) the PDB component library,
        augmented with the CYM/LYN variants pDynamo3 doesn't ship (see
        module docstring).
    """
    global _LIBRARY
    if _LIBRARY is None:
        from pBabel import PDBComponentLibrary, PDBComponentVariant
        library = PDBComponentLibrary.WithOptions()
        cym = PDBComponentVariant.WithOptions(
            componentLabel='CYS', label='Deprotonated',
            atomsToDelete=['HG'], formalCharges={'SG': -1})
        lyn = PDBComponentVariant.WithOptions(
            componentLabel='LYS', label='Neutral',
            atomsToDelete=['HZ3'], formalCharges={'NZ': 0})
        library.variants[cym.key] = cym
        library.variants[lyn.key] = lyn
        _LIBRARY = library
    return _LIBRARY


def _template_atoms_and_bonds(component_label, variant_label):
    """ Returns (atoms, bonds, applied_variant_label) for `component_label`
        with `variant_label` applied (None -> the component's own default,
        e.g. HIS -> "Delta Protonated", or the plain base template if it
        has no default, e.g. CYS/LYS/ALA/...).

        atoms: {atom_name: atomicNumber}. bonds: [(name1, name2), ...]
        (names only -- bond order/aromaticity don't matter for deciding
        which heavy atom a missing hydrogen attaches to, which is all
        this tool needs, so formalCharges/bondTypes from the variant are
        intentionally NOT applied here).

        Mirrors pBabel.PDBModel.PDBModelComponent.ApplyLibraryVariant's
        own algorithm (pBabel/PDBModel.py:167-196) -- bondsToDelete,
        atomsToDelete, atomsToAdd, bondsToAdd, in that order -- applied
        to plain name-keyed structures instead of the full PDBModel
        object graph, since our input is an already-loaded in-memory
        System, not a fresh PDB file read.

        Returns (None, None, None) if `component_label` has no library
        entry (not a standard residue).
    """
    library = _get_library()
    component = library.GetComponent(component_label)
    if component is None:
        return None, None, None

    atoms = {a.label: a.atomicNumber for a in component.atoms}
    bonds = [(b.atomLabel1, b.atomLabel2) for b in component.bonds]

    if variant_label is None:
        defaults = getattr(component, 'variants', None)
        variant_label = defaults[0] if defaults else None
    elif variant_label in _BASE_VARIANT_LABELS:
        variant_label = None

    applied = None
    if variant_label is not None:
        variant = library.GetVariant(variant_label, component_label)
        if variant is not None:
            applied = variant_label
            for pair in (variant.bondsToDelete or []):
                l1, l2 = pair
                bonds = [b for b in bonds if set(b) != {l1, l2}]
            for label in (variant.atomsToDelete or []):
                atoms.pop(label, None)
                bonds = [b for b in bonds if label not in b]
            for atom in (variant.atomsToAdd or []):
                atoms[atom.label] = atom.atomicNumber
            for bond in (variant.bondsToAdd or []):
                bonds.append((bond.atomLabel1, bond.atomLabel2))

    for name in _TERMINUS_ONLY_ATOM_NAMES:
        atoms.pop(name, None)
        bonds = [b for b in bonds if name not in b]

    return atoms, bonds, applied


def analyze_residue(component_label, existing_atom_names, variant_label=None, existing_h_bond_counts=None):
    """ Diffs the real template (see _template_atoms_and_bonds) for
        `component_label`/`variant_label` against `existing_atom_names`.

        Returns None if `component_label` isn't a standard residue
        (caller should report it as skipped), else:
          {
            'applied_variant':      str or None,
            'missing_hydrogens':    [(name, atomicNumber, parent_name), ...],
            'unbuildable':          [name, ...],  # template atoms not
                                                   # present and not
                                                   # addable (missing
                                                   # heavy atoms, or a
                                                   # missing H whose own
                                                   # parent heavy atom is
                                                   # ALSO missing)
          }

        `existing_h_bond_counts`, if given, is {parent_name: n} -- the
        number of hydrogens ALREADY REALLY BONDED (by actual geometry,
        not name) to each heavy atom in the file. Required for a
        correct result whenever the file's own hydrogen NAMES don't
        match pDynamo3's own naming convention (confirmed a real,
        common case: a file using AMBER-style "HB1/HB2" for a CH2,
        where pDynamo3's own template names the same physical pair
        "HB2/HB3" -- name-only diffing saw "HB3" absent and added a
        THIRD hydrogen onto an already-full carbon, which
        BuildHydrogenCoordinates3FromConnectivity then correctly refused
        to place, since its coordination check does look at the real,
        already-correct bond count). Without it (None), this falls back
        to counting only EXACT name matches among the template's own H
        names for that parent -- correct only when the file already
        uses pDynamo3's own naming.
    """
    atoms, bonds, applied = _template_atoms_and_bonds(component_label, variant_label)
    if atoms is None:
        return None
    existing = set(existing_atom_names)

    expected_by_parent = {}
    for name, atomic_number in atoms.items():
        if atomic_number != 1:
            continue
        parent = None
        for (l1, l2) in bonds:
            if l1 == name and l2 != name:
                parent = l2
                break
            if l2 == name and l1 != name:
                parent = l1
                break
        expected_by_parent.setdefault(parent, []).append(name)

    missing_hydrogens = []
    unbuildable = []
    for name, atomic_number in atoms.items():
        if atomic_number != 1 and name not in existing:
            unbuildable.append(name)

    for parent, h_names in expected_by_parent.items():
        if parent is None or parent not in existing:
            unbuildable.extend(n for n in h_names if n not in existing)
            continue
        if existing_h_bond_counts is not None:
            already_have = existing_h_bond_counts.get(parent, 0)
        else:
            already_have = sum(1 for n in h_names if n in existing)
        deficit = max(0, len(h_names) - already_have)
        candidates = [n for n in h_names if n not in existing]
        for name in candidates[:deficit]:
            missing_hydrogens.append((name, 1, parent))

    return {
        'applied_variant': applied,
        'missing_hydrogens': missing_hydrogens,
        'unbuildable': unbuildable,
    }


def detect_ambiguity(component_label, existing_atom_names):
    """ Returns None if there's no real protonation choice to make for
        this residue (no listed choices, or the atoms already on file
        unambiguously match exactly one), else the list of variant
        labels (from PROTONATION_VARIANTS) that are ALL still
        consistent with what's on file -- genuinely ambiguous, per the
        rule "never silently choose when there's a real choice, but
        don't ask when the file already disambiguates it."

        A variant is "consistent" iff applying it never requires
        DELETING an atom that is already present on file (e.g. HIS with
        HE2 present rules out "Delta Protonated", which deletes HE2).
    """
    choices = PROTONATION_VARIANTS.get(component_label)
    if not choices:
        return None
    existing = set(existing_atom_names)
    base_atoms, _, _ = _template_atoms_and_bonds(component_label, None if 'Protonated' not in choices else 'Protonated')
    # . Use the FULL (no-deletions) atom set as the reference for "was
    # this atom deleted by variant X" -- the union of every choice's own
    # atoms is safest, since different variants can each delete
    # different atoms relative to each other, not just relative to one
    # arbitrarily-picked "base".
    full_atoms = set()
    per_choice_atoms = {}
    for label in choices:
        variant_label = None if label in _BASE_VARIANT_LABELS else label
        atoms, _, _ = _template_atoms_and_bonds(component_label, variant_label)
        if atoms is None:
            continue
        per_choice_atoms[label] = set(atoms.keys())
        full_atoms.update(atoms.keys())
    consistent = []
    for label, atom_names in per_choice_atoms.items():
        deleted_relative_to_full = full_atoms - atom_names
        if existing & deleted_relative_to_full:
            continue  # this choice would have deleted an atom that IS present -> contradicted
        consistent.append(label)
    if len(consistent) <= 1:
        return None
    return consistent


def build_residue_index(system):
    """ Groups `system`'s atoms into residues, in encounter order.

        Returns (residues, atom_index_to_key):
          residues: [ (key, [atom_name, ...]), ... ] in encounter order,
                     key = (chain, resName, resSeq)
          atom_index_to_key: {pdynamo_atom.index: key}

        Uses the same atom.parent/atom.parent.parent + Sequence.ParseLabel
        pattern EasyHybrid's own
        _get_atom_info_from_pdynamo_atom_obj()
        (pdynamo/pDynamo2EasyHybrid/session.py) already relies on.
    """
    sequence = getattr(system, 'sequence', None)
    residues = {}
    order = []
    atom_index_to_key = {}
    for atom in system.atoms:
        entity_label = atom.parent.parent.label
        chain = entity_label[0:1]
        if sequence is not None:
            resName, resSeq, _iCode = sequence.ParseLabel(atom.parent.label, fields=3)
        else:
            resName, resSeq = atom.parent.label, '1'
        key = (chain, resName, resSeq)
        if key not in residues:
            residues[key] = []
            order.append(key)
        residues[key].append(atom.label)
        atom_index_to_key[atom.index] = key
    return [(key, residues[key]) for key in order], atom_index_to_key


def build_residue_atom_objects(system):
    """ Same grouping as build_residue_index(), but keyed to the actual
        Atom OBJECTS of each residue (needed for real-connectivity-based
        checks, e.g. _existing_h_bond_counts_for_residue()) instead of
        just their names.

        Returns {key: [Atom, ...]}, key = (chain, resName, resSeq).
    """
    sequence = getattr(system, 'sequence', None)
    residues_map = {}
    for atom in system.atoms:
        entity_label = atom.parent.parent.label
        chain = entity_label[0:1]
        if sequence is not None:
            resName, resSeq, _iCode = sequence.ParseLabel(atom.parent.label, fields=3)
        else:
            resName, resSeq = atom.parent.label, '1'
        key = (chain, resName, resSeq)
        residues_map.setdefault(key, []).append(atom)
    return residues_map


def _ensure_connectivity(system):
    """ Makes sure `system.connectivity` reflects its REAL geometry --
        a freshly-imported PDB with no mmModel defined yet may have no
        connectivity at all, and even one with mmState-derived bonds is
        fine to use as-is (BondsFromCoordinates3() only acts when
        connectivity is missing/empty).
    """
    if system.connectivity is None or len(system.connectivity.bonds) == 0:
        system.BondsFromCoordinates3()


def _existing_h_bond_counts_for_residue(system, atoms_in_residue):
    """ {heavy_atom_label: n} -- for one residue's real Atom objects,
        how many hydrogen neighbors each heavy atom ALREADY has, per
        `system.connectivity`'s real (geometry-derived) bonds -- NOT
        atom names. See analyze_residue()'s own docstring for why this
        matters: a file using a different H-naming convention than
        pDynamo3's own template still has the CORRECT number of real
        hydrogen bonds, and this is what must be trusted, not names.
    """
    adjacent = system.connectivity.adjacentNodes
    counts = {}
    for atom in atoms_in_residue:
        if atom.atomicNumber == 1:
            continue
        counts[atom.label] = sum(1 for neighbor in adjacent.get(atom, ()) if neighbor.atomicNumber == 1)
    return counts


# . Generous upper bound (Angstrom) for a real disulfide S-S bond
# (~2.05 A) -- wide enough to tolerate a slightly strained/modelled
# bridge, far short of any plausible non-bonded S...S contact.
_DISULFIDE_SG_DISTANCE = 2.5


def find_disulfide_bonded_cysteines(system):
    """ Returns the set of (chain, resName, resSeq) keys for every CYS
        residue whose own SG atom sits within real disulfide-bond
        distance of another CYS's SG.

        These residues are chemically bonded, not free thiols -- they
        must never get a new HG (nor be offered a protonated/
        deprotonated choice, which only makes sense for a free -SH/-S-).
        pDynamo3's own component library models this as a distinct
        "Disulfide Bridge" LINK (pBabel/PDBComponentLibrary.py:169-170)
        that deletes HG from BOTH partners when reading a fresh PDB file
        with explicit links -- this reproduces that same exclusion by
        plain geometry instead, since the caller here is an
        already-loaded, in-memory System, not a fresh PDB-with-links
        read. Found necessary in practice: a real test fixture
        (crambin.pdb) has all 6 of its cysteines disulfide-bonded (3
        bridges), and without this check every one of them would be
        flagged as an "ambiguous protonated/deprotonated" free thiol.
    """
    sg_atoms = []
    sequence = getattr(system, 'sequence', None)
    for atom in system.atoms:
        if atom.label != 'SG':
            continue
        entity_label = atom.parent.parent.label
        chain = entity_label[0:1]
        if sequence is not None:
            resName, resSeq, _iCode = sequence.ParseLabel(atom.parent.label, fields=3)
        else:
            resName, resSeq = atom.parent.label, '1'
        if resName == 'CYS':
            sg_atoms.append(((chain, resName, resSeq), atom))
    bonded = set()
    coords = system.coordinates3
    cutoff2 = _DISULFIDE_SG_DISTANCE * _DISULFIDE_SG_DISTANCE
    for i in range(len(sg_atoms)):
        key_i, atom_i = sg_atoms[i]
        xi, yi, zi = coords[atom_i.index, 0], coords[atom_i.index, 1], coords[atom_i.index, 2]
        for j in range(i + 1, len(sg_atoms)):
            key_j, atom_j = sg_atoms[j]
            xj, yj, zj = coords[atom_j.index, 0], coords[atom_j.index, 1], coords[atom_j.index, 2]
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            if d2 <= cutoff2:
                bonded.add(key_i)
                bonded.add(key_j)
    return bonded


def analyze_system(system, residue_keys=None):
    """ Analyzes every residue in `system` (or only those in
        `residue_keys`, for the "selected atoms" scope) and returns a
        list of per-residue result dicts -- only for residues that are
        either non-standard (no library entry) or have 1+ missing
        hydrogens/unbuildable atoms. Residues needing nothing are
        omitted entirely.

        Each dict:
          {
            'key': (chain, resName, resSeq),
            'component_label': resName,
            'is_standard_residue': bool,
            'ambiguous_choices': [str, ...] or None,
            'default_choice': str or None,   # 'Protonated'/'Deprotonated'/...
            'missing_hydrogens': [(name, atomicNumber, parent_name), ...],
            'unbuildable': [name, ...],
            'disulfide_bonded': bool,   # CYS only -- see find_disulfide_bonded_cysteines()
          }
    """
    _ensure_connectivity(system)
    residues, _ = build_residue_index(system)
    residue_atoms = build_residue_atom_objects(system)
    if residue_keys is not None:
        wanted = set(residue_keys)
        residues = [(k, names) for (k, names) in residues if k in wanted]
    disulfide_keys = find_disulfide_bonded_cysteines(system)
    results = []
    for (key, names) in residues:
        chain, resName, resSeq = key
        is_disulfide = key in disulfide_keys
        h_counts = _existing_h_bond_counts_for_residue(system, residue_atoms.get(key, []))
        analysis = analyze_residue(resName, names, variant_label=None, existing_h_bond_counts=h_counts)
        if analysis is None:
            # . No PDBComponentLibrary entry -- a ligand/cofactor/any
            # other organic molecule. Not run through OpenBabel here
            # (that's an external-process round trip PER residue --
            # deferred to commit time in
            # rebuild_system_with_added_hydrogens(), so analyzing a
            # whole system doesn't pay that cost for residues the user
            # might not even keep in scope). Just report enough for the
            # GUI to offer a pH field instead of silently skipping.
            n_heavy = sum(1 for a in residue_atoms.get(key, []) if a.atomicNumber != 1)
            results.append({
                'key': key, 'component_label': resName, 'is_standard_residue': False,
                'is_organic_candidate': n_heavy > 0,
                'n_heavy_atoms': n_heavy,
            })
            continue
        missing_hydrogens = analysis['missing_hydrogens']
        if is_disulfide:
            # . A disulfide-bonded SG is never a free thiol -- HG must
            # never be added regardless of any protonated/deprotonated
            # choice.
            missing_hydrogens = [m for m in missing_hydrogens if m[0] != 'HG']
        if not missing_hydrogens and not analysis['unbuildable']:
            continue
        ambiguous = None if is_disulfide else detect_ambiguity(resName, names)
        default_choice = ambiguous[0] if ambiguous else (analysis['applied_variant'] or 'Protonated')
        results.append({
            'key': key,
            'component_label': resName,
            'is_standard_residue': True,
            'disulfide_bonded': is_disulfide,
            'ambiguous_choices': ambiguous,
            'default_choice': default_choice,
            'missing_hydrogens': missing_hydrogens,
            'unbuildable': analysis['unbuildable'],
        })
    return results


def rebuild_system_with_added_hydrogens(old_system, residue_variant_choices, random_seed=None):
    """ Builds a NEW pDynamo System containing every atom `old_system`
        already had (same Atom objects, same old coordinates) PLUS the
        missing hydrogens for each residue named in
        `residue_variant_choices` (a {(chain, resName, resSeq):
        variant_label_or_None} mapping -- typically built from
        analyze_system()'s own 'key'/'default_choice' fields, after the
        caller lets the user review/edit ambiguous ones).

        New hydrogen coordinates are built by
        pSimulation.BuildHydrogenCoordinates3FromConnectivity, using only
        the (old + newly-declared) bond graph -- no forcefield involved.

        Returns a report dict:
          {
            'system': new System, or None if nothing needed adding,
            'added': {key: [new_atom_name, ...]},
            'skipped_non_standard': [key, ...],
            'unbuildable': {key: [name, ...]},
            'unbuilt_by_geometry': [name, ...],  # still undefined after
                                                  # BuildHydrogenCoordinates3FromConnectivity
                                                  # -- never silently
                                                  # dropped, always
                                                  # reported by name.
          }
    """
    from pMolecule import Atom, System, BondType
    from pMolecule.Sequence import Sequence
    from pScientific.Geometry3 import Coordinates3
    from pScientific.RandomNumbers import RandomNumberGenerator
    from pSimulation import BuildHydrogenCoordinates3FromConnectivity

    _ensure_connectivity(old_system)
    residues, _ = build_residue_index(old_system)
    residue_order = [key for key, _names in residues]

    # . Group actual Atom OBJECTS (not just names) per residue key, in
    # the same order build_residue_index used, so downstream ordering
    # of atom_paths/atoms_list stays consistent with `residues`.
    residues_map = build_residue_atom_objects(old_system)

    report = {'added': {}, 'skipped_non_standard': [], 'unbuildable': {}, 'charge_changes': {}}
    disulfide_keys = find_disulfide_bonded_cysteines(old_system)
    residue_additions = {}
    atoms_to_remove = set()
    for key, variant_choice in residue_variant_choices.items():
        atoms_in_residue = residues_map.get(key)
        if not atoms_in_residue:
            continue

        if isinstance(variant_choice, tuple) and variant_choice[0] == ORGANIC_CHOICE_KIND:
            ph = variant_choice[1]
            organic_result = protonate_non_standard_residue(old_system, atoms_in_residue, ph=ph)
            if organic_result is None:
                report['skipped_non_standard'].append(key)
                continue
            additions = [(name, 1, parent_atom, xyz)
                         for (name, _an, parent_atom, xyz) in organic_result['missing_hydrogens']]
            if organic_result['charge_changes']:
                report['charge_changes'][key] = organic_result['charge_changes']
            # . Any hydrogens the residue already had are replaced
            # outright by OpenBabel's own fresh, complete set above --
            # see protonate_non_standard_residue()'s own docstring for
            # why keeping both would silently double-count valence.
            atoms_to_remove.update(id(a) for a in organic_result['atoms_to_remove'])
            if additions:
                residue_additions[key] = additions
                report['added'][key] = [a[0] for a in additions]
            continue

        chain, resName, resSeq = key
        existing_names = [a.label for a in atoms_in_residue]
        h_counts = _existing_h_bond_counts_for_residue(old_system, atoms_in_residue)
        analysis = analyze_residue(resName, existing_names, variant_label=variant_choice, existing_h_bond_counts=h_counts)
        if analysis is None:
            report['skipped_non_standard'].append(key)
            continue
        name_to_atom = {a.label: a for a in atoms_in_residue}
        additions = []
        for (new_name, atomic_number, parent_name) in analysis['missing_hydrogens']:
            if key in disulfide_keys and new_name == 'HG':
                # . Never add a thiol H to a disulfide-bonded SG -- see
                # find_disulfide_bonded_cysteines()'s own docstring.
                continue
            parent_atom = name_to_atom.get(parent_name)
            if parent_atom is None:
                continue
            additions.append((new_name, atomic_number, parent_atom, None))
        if additions:
            residue_additions[key] = additions
            report['added'][key] = [a[0] for a in additions]
        if analysis['unbuildable']:
            report['unbuildable'][key] = list(analysis['unbuildable'])

    if not residue_additions:
        report['system'] = None
        report['unbuilt_by_geometry'] = []
        return report

    # . old_system.connectivity was already ensured at the top of this
    # function (needed for the real-bond-count check above too) -- carry
    # it forward as-is into the new system's own bond list below.

    # . Capture every OLD atom's real (x,y,z) BEFORE anything reassigns
    # .index -- Sequence.FromAtomPaths()/System.FromSequence() both
    # reassign every atom's .index (see vismol_pdbqt_builder.py's own
    # documented fix for this exact footgun: writing coordinates by
    # ORIGINAL index instead of the atom's FINAL index silently swaps
    # atoms' displayed positions).
    old_coords = old_system.coordinates3
    captured_xyz = {}
    for atom in old_system.atoms:
        i = atom.index
        captured_xyz[id(atom)] = (old_coords[i, 0], old_coords[i, 1], old_coords[i, 2])

    # . Build the new full atom list, REUSING the same Atom objects for
    # everything that already existed -- so old_system.connectivity's own
    # Bond objects (which reference those same Atom objects) stay valid
    # with no translation needed, and only the newly-added H atoms need
    # brand-new Atom objects + bonds.
    atoms_list = []
    atom_paths = []
    xyz_list = []
    new_bonds = []
    for key in residue_order:
        chain, resName, resSeq = key
        for atom in residues_map.get(key, []):
            if id(atom) in atoms_to_remove:
                continue
            atoms_list.append(atom)
            atom_paths.append('{}:{}.{}:{}'.format(chain, resName, resSeq, atom.label))
            xyz_list.append(captured_xyz[id(atom)])
        for (new_name, atomic_number, parent_atom, xyz) in residue_additions.get(key, []):
            new_atom = Atom.WithOptions(atomicNumber=atomic_number, label=new_name)
            atoms_list.append(new_atom)
            atom_paths.append('{}:{}.{}:{}'.format(chain, resName, resSeq, new_name))
            xyz_list.append(xyz)
            new_bonds.append((new_atom, parent_atom, BondType.Single))

    sequence = Sequence.FromAtomPaths(atom_paths, atoms=atoms_list)
    # . Drop any old bond touching a removed atom (a pre-existing H on a
    # re-protonated organic residue, see atoms_to_remove above) --
    # Connectivity/Bond.FromIterable requires both endpoints to be
    # present in the new atom list.
    old_bonds = [b for b in old_system.connectivity.bonds
                 if id(b.node1) not in atoms_to_remove and id(b.node2) not in atoms_to_remove]
    all_bonds = old_bonds + new_bonds
    new_system = System.FromSequence(sequence, bonds=all_bonds)
    new_system.label = getattr(old_system, 'label', None)

    coords = Coordinates3.WithExtent(len(atoms_list))
    for i, atom in enumerate(atoms_list):
        row = atom.index
        xyz = xyz_list[i]
        if xyz is None:
            coords.FlagCoordinateAsUndefined(row)
        else:
            coords[row, 0], coords[row, 1], coords[row, 2] = xyz
    new_system.coordinates3 = coords

    if random_seed is not None:
        rng = RandomNumberGenerator.WithSeed(random_seed)
    else:
        rng = RandomNumberGenerator.WithRandomSeed()
    BuildHydrogenCoordinates3FromConnectivity(new_system, randomNumberGenerator=rng)

    unbuilt = []
    if new_system.coordinates3.numberUndefined > 0:
        index_to_atom = {a.index: a for a in new_system.atoms}
        undefined_rows = getattr(new_system.coordinates3, 'undefined', None) or []
        for row in sorted(undefined_rows):
            a = index_to_atom.get(row)
            if a is not None:
                unbuilt.append(a.label)

    report['system'] = new_system
    report['unbuilt_by_geometry'] = unbuilt
    return report


# =============================================================================
#  Non-standard residues / ligands / organic molecules -- OpenBabel-based
# =============================================================================
# . A residue name with no PDBComponentLibrary entry has no template to
# diff against -- there is no "expected atom list" to look up. OpenBabel
# (already a hard dependency of this app, already used identically for
# ligand docking prep -- vina_runner.convert_ligand_file_to_pdbqt_ensemble()'s
# own `-p <ph>` flag) solves this for an arbitrary organic molecule from
# geometry alone: it perceives bonds/valence from the heavy-atom
# positions and adds hydrogens with a real, pH-dependent ionization
# state via its own per-atom pKa model -- verified directly this
# session against a real ligand fixture (DB00220.pdb from this
# project's own docking test data): stripping its hydrogens and
# re-running `obabel ... -p <ph>` reproduced the original heavy-atom
# positions/order EXACTLY (0 mismatches), and comparing `-h` (plain
# valence) against `-p 7.4` showed a real extra proton on a basic
# nitrogen at physiological pH -- confirming this is genuine
# pH-dependent protonation modeling, not just valence-filling.

# . Marks a residue_variant_choices value as "use OpenBabel, not
# PDBComponentLibrary" -- ('organic', ph) instead of a variant-label
# string/None.
ORGANIC_CHOICE_KIND = 'organic'

_PDB_ATOM_LINE = "{:<6s}{:5d} {:<4s} {:3s} {:1s}{:>4s}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
# ^ same fixed-column format as pdynamo/pDynamo2EasyHybrid/helpers.py's
# own export_special_PDB() -- already fixed once this session for a
# real OpenBabel-strictness bug (columns 77-78 must be the RIGHT-
# justified 2-char element symbol), reused verbatim rather than
# re-deriving it.


def _write_minimal_pdb(atoms_with_xyz, path, resname='LIG', chain='A', resseq=1):
    """ Writes a minimal, single-residue PDB (ATOM records only) for
        `atoms_with_xyz` -- a list of (element_symbol, x, y, z).
        resname/chain/resseq are cosmetic placeholders for OpenBabel's
        own benefit (it only uses element + geometry to perceive bonds/
        valence) -- atom NAMES are deliberately NOT taken from the real
        system here, since matching the result back to the real atoms
        is done by POSITION (confirmed reliable -- see module docstring
        above), not by name, so a placeholder "X1"/"X2"/... naming
        scheme is used for the temp file only.
    """
    with open(path, 'w') as f:
        for i, (element, x, y, z) in enumerate(atoms_with_xyz, start=1):
            name = '{}{}'.format(element, i)
            f.write(_PDB_ATOM_LINE.format(
                'ATOM  ', i, name, resname, chain, str(resseq), x, y, z, 1.0, 0.0, element))
        f.write('END\n')


def _parse_pdb_atoms_and_conect(path):
    """ Parses a plain PDB written/read by OpenBabel -- returns
        (atoms, conect) where atoms is an ordered list of
        {serial, element, formal_charge, x, y, z} (1-based serial,
        matching input order) and conect is {serial: [other_serial, ...]}
        from every CONECT record (both directions, as OpenBabel writes
        them redundantly -- confirmed in a real test file).
    """
    atoms = []
    conect = {}
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')):
            serial = int(line[6:11])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip()
            charge_str = line[78:80].strip()
            formal_charge = 0
            if charge_str:
                digit, sign = charge_str[0], charge_str[1:2]
                if digit.isdigit():
                    formal_charge = int(digit) * (-1 if sign == '-' else 1)
            atoms.append({'serial': serial, 'element': element,
                           'formal_charge': formal_charge, 'x': x, 'y': y, 'z': z})
        elif line.startswith('CONECT'):
            parts = line.split()
            if len(parts) >= 2:
                serial = int(parts[1])
                others = [int(p) for p in parts[2:]]
                conect.setdefault(serial, []).extend(others)
    return atoms, conect


def protonate_non_standard_residue(system, atoms_in_residue, ph=7.4, obabel_bin=None):
    """ Adds missing hydrogens AND assigns pH-dependent protonation
        state to a residue with NO PDBComponentLibrary entry (a ligand,
        cofactor, or any other organic molecule), via OpenBabel -- see
        this section's own header comment for why/how this differs from
        the standard-residue path.

        `atoms_in_residue`: real pDynamo Atom objects (same ones
        build_residue_atom_objects()/build_residue_index() group by
        residue key) -- their CURRENT heavy-atom geometry is exported,
        OpenBabel adds/ionizes hydrogens, and the result is mapped back
        by POSITION (heavy atoms are 1:1, in the same order, confirmed
        stable) plus each new hydrogen's own CONECT record (maps it to
        its parent's serial number, <= the heavy atom count).

        Returns None if OpenBabel isn't available, or a dict shaped
        like analyze_residue()'s own:
          {
            'missing_hydrogens': [(new_name, 1, parent_atom_object, (x,y,z)), ...],
            'charge_changes': [(existing_atom_label, new_formal_charge), ...],
            'atoms_to_remove': [atom_object, ...],
          }
        `missing_hydrogens` here carries the PARENT ATOM OBJECT directly
        (not a name needing a further lookup) and an explicit (x,y,z) --
        both differences from analyze_residue()'s own shape, needed
        because there's no template/name-index to resolve against, and
        because OpenBabel already computed real 3D positions (no
        BuildHydrogenCoordinates3FromConnectivity call needed for these).

        Only HEAVY atoms currently on the residue are sent to OpenBabel
        -- any hydrogens already present are excluded from the INPUT and
        `atoms_to_remove` lists them so the caller
        (rebuild_system_with_added_hydrogens()) drops them from the
        final atom list too, letting OpenBabel's own fresh, complete set
        (in `missing_hydrogens`) replace them outright -- avoiding a
        mixed "some H's ours (kept), some pre-existing (also kept),
        some new (added)" double-count that would silently over-value
        every protonated heavy atom's valence.

        New hydrogen names are guaranteed <= 4 characters (the fixed
        PDB atom-name field width this codebase's own PDB writers use,
        e.g. helpers.py's export_special_PDB ATOMLINE format,
        `"{:<4s}"` -- a longer string isn't truncated, it just breaks
        column alignment in any later export). Preferred form is
        "H"+parent's own label when that already fits (e.g. parent "C1"
        -> "HC1"); anything that wouldn't fit (a long parent label, or a
        parent needing more than one new H) falls back to a plain
        "H"+running-index scheme, since these residues have no
        canonical per-atom naming convention to preserve anyway.
    """
    from pScientific import PeriodicTable

    obabel_bin = obabel_bin or shutil.which('obabel')
    if not obabel_bin:
        return None

    heavy_atoms = [a for a in atoms_in_residue if a.atomicNumber != 1]
    existing_h_atoms = [a for a in atoms_in_residue if a.atomicNumber == 1]
    if not heavy_atoms:
        return None

    atoms_with_xyz = []
    for atom in heavy_atoms:
        i = atom.index
        coords = system.coordinates3
        element = PeriodicTable.Symbol(atom.atomicNumber)
        atoms_with_xyz.append((element, coords[i, 0], coords[i, 1], coords[i, 2]))

    with tempfile.TemporaryDirectory(prefix='eh3_protonate_') as tmp_dir:
        in_path = os.path.join(tmp_dir, 'in.pdb')
        out_path = os.path.join(tmp_dir, 'out.pdb')
        _write_minimal_pdb(atoms_with_xyz, in_path)
        cmd = [obabel_bin, in_path, '-O', out_path, '-p', str(ph)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not os.path.exists(out_path):
            return None
        out_atoms, conect = _parse_pdb_atoms_and_conect(out_path)

    n_heavy = len(heavy_atoms)
    if len(out_atoms) < n_heavy:
        # . OpenBabel produced fewer atoms than it was given -- something
        # went chemically wrong (e.g. it couldn't parse the geometry at
        # all); refuse rather than guess at a broken correspondence.
        return None

    charge_changes = []
    for i, atom in enumerate(heavy_atoms):
        new_charge = out_atoms[i]['formal_charge']
        if new_charge != getattr(atom, 'formalCharge', 0):
            atom.formalCharge = new_charge
            charge_changes.append((atom.label, new_charge))

    new_h_by_parent = {}
    for out_atom in out_atoms[n_heavy:]:
        if out_atom['element'] != 'H':
            continue
        parents = [p for p in conect.get(out_atom['serial'], []) if p <= n_heavy]
        if not parents:
            continue
        parent_atom = heavy_atoms[parents[0] - 1]
        new_h_by_parent.setdefault(id(parent_atom), []).append((parent_atom, out_atom))

    # . Name new H's -- ALWAYS <= 4 characters (see this function's own
    # docstring for why this is a hard requirement, not cosmetic).
    # Preferred: "H" + parent's own label, when a parent gets exactly
    # one new H AND that name is short enough to fit. Fallback: a plain
    # "H"+running-index, unique within the residue -- used whenever the
    # preferred form would exceed 4 characters (a long parent label, or
    # more than one new H on the same parent, which would otherwise need
    # a "_1"/"_2" suffix that rarely fits).
    existing_names = {a.label for a in atoms_in_residue}
    used_names = set(existing_names)
    fallback_counter = [1]

    def _next_fallback_name():
        while True:
            candidate = 'H{}'.format(fallback_counter[0])
            fallback_counter[0] += 1
            if len(candidate) <= 4 and candidate not in used_names:
                return candidate

    missing_hydrogens = []
    for entries in new_h_by_parent.values():
        multiple = len(entries) > 1
        for (parent_atom, out_atom) in entries:
            name = None
            if not multiple:
                candidate = 'H' + parent_atom.label
                if len(candidate) <= 4 and candidate not in used_names:
                    name = candidate
            if name is None:
                name = _next_fallback_name()
            used_names.add(name)
            xyz = (out_atom['x'], out_atom['y'], out_atom['z'])
            missing_hydrogens.append((name, 1, parent_atom, xyz))

    return {
        'missing_hydrogens': missing_hydrogens,
        'charge_changes': charge_changes,
        'atoms_to_remove': existing_h_atoms,
    }


# =============================================================================
#  Remove Hydrogens
# =============================================================================
# . The complementary operation to everything above -- no residue
# template, no OpenBabel, no geometry building needed: for each real
# hydrogen atom already in the system, decide (by its own real bonded
# neighbor, from actual connectivity, not by name) whether to remove it,
# then rebuild the System excluding the chosen ones. Same fixed-size-
# container constraint as the addition path (System.atoms/coordinates3
# can't be resized in place), so this still needs a full
# Sequence.FromAtomPaths + System.FromSequence rebuild -- just a much
# simpler one, with nothing to ADD.

# . Apolar/nonpolar = bonded to carbon (C-H bonds -- C and H have
# similar electronegativity, no meaningful bond dipole). Anything else
# (N-H, O-H, S-H, ...) is polar -- it can participate in hydrogen
# bonding and is kept by the 'apolar' mode. This is the same
# "united-atom" convention AutoDock/AutoDockTools' own "merge nonpolar
# hydrogens" step uses (this app's own docking tools already rely on
# the same distinction conceptually, via OpenBabel's `--AddPolarH` for
# receptors).
_APOLAR_PARENT_ATOMIC_NUMBER = 6  # carbon


def find_removable_hydrogens(system, mode='all', residue_keys=None):
    """ Returns the list of real hydrogen Atom objects in `system`
        matching `mode`:
          'all'    -- every hydrogen atom in scope.
          'apolar' -- only hydrogens bonded to a carbon (see this
                      section's own header comment).

        `residue_keys`, if given, restricts scope to those residues
        only (the "selected atoms" scope, same residue-level
        granularity build_residue_index()/analyze_system() already use
        elsewhere in this module) -- None means the entire system.
    """
    _ensure_connectivity(system)
    _residues, atom_index_to_key = build_residue_index(system)
    wanted = set(residue_keys) if residue_keys is not None else None
    adjacent = system.connectivity.adjacentNodes
    to_remove = []
    for atom in system.atoms:
        if atom.atomicNumber != 1:
            continue
        if wanted is not None and atom_index_to_key.get(atom.index) not in wanted:
            continue
        if mode == 'apolar':
            neighbors = adjacent.get(atom, ())
            if not any(n.atomicNumber == _APOLAR_PARENT_ATOMIC_NUMBER for n in neighbors):
                continue
        to_remove.append(atom)
    return to_remove


def remove_hydrogens_from_system(old_system, atoms_to_remove):
    """ Builds a new pDynamo System excluding `atoms_to_remove` (real
        Atom objects from old_system.atoms) -- the remove-only
        counterpart to rebuild_system_with_added_hydrogens(), reusing
        the same core assembly pattern (real-coordinate capture by
        object identity, Sequence.FromAtomPaths + System.FromSequence,
        old-bond filtering) but with none of that function's
        addition-specific machinery (no new atoms, no
        BuildHydrogenCoordinates3FromConnectivity call needed).

        Returns the new System (never None -- even removing zero atoms
        would just rebuild an equivalent system; callers should check
        `atoms_to_remove` themselves before calling this, same as
        find_removable_hydrogens() returning an empty list).
    """
    from pMolecule import System
    from pMolecule.Sequence import Sequence
    from pScientific.Geometry3 import Coordinates3

    _ensure_connectivity(old_system)
    residues, _ = build_residue_index(old_system)
    residue_order = [key for key, _names in residues]
    residues_map = build_residue_atom_objects(old_system)
    to_remove_ids = {id(a) for a in atoms_to_remove}

    old_coords = old_system.coordinates3
    captured_xyz = {}
    for atom in old_system.atoms:
        i = atom.index
        captured_xyz[id(atom)] = (old_coords[i, 0], old_coords[i, 1], old_coords[i, 2])

    atoms_list = []
    atom_paths = []
    xyz_list = []
    for key in residue_order:
        chain, resName, resSeq = key
        for atom in residues_map.get(key, []):
            if id(atom) in to_remove_ids:
                continue
            atoms_list.append(atom)
            atom_paths.append('{}:{}.{}:{}'.format(chain, resName, resSeq, atom.label))
            xyz_list.append(captured_xyz[id(atom)])

    sequence = Sequence.FromAtomPaths(atom_paths, atoms=atoms_list)
    old_bonds = [b for b in old_system.connectivity.bonds
                 if id(b.node1) not in to_remove_ids and id(b.node2) not in to_remove_ids]
    new_system = System.FromSequence(sequence, bonds=old_bonds)
    new_system.label = getattr(old_system, 'label', None)

    coords = Coordinates3.WithExtent(len(atoms_list))
    for i, atom in enumerate(atoms_list):
        row = atom.index
        x, y, z = xyz_list[i]
        coords[row, 0], coords[row, 1], coords[row, 2] = x, y, z
    new_system.coordinates3 = coords

    return new_system
