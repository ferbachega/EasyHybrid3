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
    "OXT" and "HXT" are always excluded from an amino acid's expected atom
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

from util.chain_ids import system_chain_id_map, chain_id_for_atom

# . Atom names that only ever belong to a genuine chain terminus -- see
# the module docstring's SCOPE note. Never added by this tool.
_TERMINUS_ONLY_ATOM_NAMES = frozenset(('H2', 'OXT', 'HXT'))

# [EN] 2026-10-02 user request ("para metais, o padrao e nao adicionar
# hidrogenios"): metal atoms never receive hydrogens by default -- neither
# here (OpenBabel would turn a Mg2+ ion into MgH2) nor in the Builder's own
# atom_ops.adjust_hydrogens(). The set is vismol's own METAL_ELEMENTS
# (alkali, alkaline-earth, transition, lanthanides, actinides -- the one that
# also drives the dashed metal-bond rendering, left untouched) plus the
# post-transition metals, which also form ions/complexes rather than hydrides.
_POST_TRANSITION_METALS = frozenset(('Al', 'Ga', 'In', 'Sn', 'Tl', 'Pb', 'Bi', 'Po'))


def is_hydrogen_free_metal(symbol):
    """ True if `symbol` (any case, e.g. 'MG', 'Zn') is a metal that should
        get no hydrogens by default (see the note above). """
    s = ''.join(ch for ch in str(symbol) if ch.isalpha())
    if not s:
        return False
    s = s[0].upper() + s[1:].lower()
    if s in _POST_TRANSITION_METALS:
        return True
    try:
        from vismol.core.metal_bonds import METAL_ELEMENTS
    except Exception:
        return False
    return s in METAL_ELEMENTS


def _atomic_number_is_hydrogen_free_metal(atomic_number):
    try:
        from pScientific import PeriodicTable
        return is_hydrogen_free_metal(PeriodicTable.Symbol(int(atomic_number)))
    except Exception:
        return False

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


def _template_atoms_and_bonds(component_label, variant_label, strip_terminus_atoms=True):
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

    # [EN] 2026-10-02 BUG FIX (user report: Add Hydrogens failed on 1BX4.pdb
    # with "'NoneType' object is not iterable"): single-atom components of
    # the pDynamo library (ions such as CL, NA) have bonds = None.
    atoms = {a.label: a.atomicNumber for a in (component.atoms or [])}
    bonds = [(b.atomLabel1, b.atomLabel2) for b in (component.bonds or [])]

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

    # [EN] 2026-10-02 BUG FIX (user report: free oxygens -- waters -- got only
    # ONE hydrogen): H2/OXT/HXT are chain-terminus-only names for AMINO ACIDS
    # (H2 = 2nd N-terminal H), but other components use the same names for
    # ordinary atoms -- HOH's second hydrogen is "H2". Only strip them from
    # amino-acid templates (those with the N/CA/C backbone).
    if strip_terminus_atoms and {'N', 'CA', 'C'} <= set(atoms):
        for name in _TERMINUS_ONLY_ATOM_NAMES:
            atoms.pop(name, None)
            bonds = [b for b in bonds if name not in b]

    return atoms, bonds, applied


def _template_capacity(component_label, variant_label):
    """ {atom_name: number of neighbours in the FULL template} -- the free
        component, terminus-only atoms (H2/OXT/HXT) included. This is the
        atom's real number of bonding partners (heavy + H), used to cap
        the hydrogens a heavy atom may receive once its bonds to OTHER
        residues are counted (see analyze_residue()). """
    atoms, bonds, _ = _template_atoms_and_bonds(component_label, variant_label, strip_terminus_atoms=False)
    if atoms is None:
        return None
    capacity = {name: 0 for name in atoms}
    for (l1, l2) in bonds:
        if l1 in capacity: capacity[l1] += 1
        if l2 in capacity: capacity[l2] += 1
    return capacity


def analyze_residue(component_label, existing_atom_names, variant_label=None, existing_h_bond_counts=None,
                    existing_heavy_bond_counts=None):
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

        `existing_heavy_bond_counts`, if given, is {parent_name: n} -- the
        number of HEAVY atoms really bonded to each heavy atom, including
        atoms of OTHER residues (peptide bonds, disulfides, links). The
        hydrogens of a parent are then capped at its full-template
        capacity minus those heavy bonds: a chain PROLINE N (bonded to CA,
        CD and the previous residue's C) gets no hydrogen -- the library's
        PRO is the free amino acid, whose N carries one H. 2026-10-04 bug
        fix: every chain proline used to receive an amide H (4-bonded N).

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
    capacity = _template_capacity(component_label, variant_label) if existing_heavy_bond_counts is not None else None

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
        if capacity is not None and parent in capacity:
            free_slots = capacity[parent] - existing_heavy_bond_counts.get(parent, 0) - already_have
            deficit = min(deficit, max(0, free_slots))
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
    chain_map = system_chain_id_map(system)
    residues = {}
    order = []
    atom_index_to_key = {}
    for atom in system.atoms:
        # 2026-09-29: same chain-ID rule as the vismol object (util/chain_ids.py);
        # entity_label[0:1] merged e.g. CHARMM segments AAAA/AABA, so same-numbered
        # residues of two different proteins were treated as ONE residue here.
        chain = chain_id_for_atom(atom, chain_map)
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
    chain_map = system_chain_id_map(system)
    residues_map = {}
    for atom in system.atoms:
        # 2026-09-29: same chain-ID rule as the vismol object (util/chain_ids.py);
        # entity_label[0:1] merged e.g. CHARMM segments AAAA/AABA, so same-numbered
        # residues of two different proteins were treated as ONE residue here.
        chain = chain_id_for_atom(atom, chain_map)
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


def _existing_heavy_bond_counts_for_residue(system, atoms_in_residue):
    """ {heavy_atom_label: n} -- number of HEAVY neighbours of each heavy
        atom of one residue, counting neighbours in other residues too
        (peptide bonds, disulfides). See analyze_residue(). """
    adjacent = system.connectivity.adjacentNodes
    counts = {}
    for atom in atoms_in_residue:
        if atom.atomicNumber == 1:
            continue
        counts[atom.label] = sum(1 for neighbor in adjacent.get(atom, ()) if neighbor.atomicNumber != 1)
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
    chain_map = system_chain_id_map(system)
    for atom in system.atoms:
        if atom.label != 'SG':
            continue
        # 2026-09-29: same chain-ID rule as the vismol object (util/chain_ids.py);
        # entity_label[0:1] merged e.g. CHARMM segments AAAA/AABA, so same-numbered
        # residues of two different proteins were treated as ONE residue here.
        chain = chain_id_for_atom(atom, chain_map)
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
        heavy_counts = _existing_heavy_bond_counts_for_residue(system, residue_atoms.get(key, []))
        analysis = analyze_residue(resName, names, variant_label=None, existing_h_bond_counts=h_counts,
                                   existing_heavy_bond_counts=heavy_counts)
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
            # metal ions (Mg2+, Zn2+, ...) get no hydrogens by default: a
            # residue made only of metal atoms is not offered to OpenBabel
            n_metal = sum(1 for a in residue_atoms.get(key, [])
                          if a.atomicNumber != 1 and _atomic_number_is_hydrogen_free_metal(a.atomicNumber))
            results.append({
                'key': key, 'component_label': resName, 'is_standard_residue': False,
                'is_organic_candidate': n_heavy - n_metal > 0,
                'is_metal_ion': n_heavy > 0 and n_metal == n_heavy,
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


# =============================================================================
#  Missing heavy atoms of amino acids (side chains) -- 2026-10-18
# =============================================================================
# . [EN] User's request: when hydrogens are added, amino-acid residues with
#   missing heavy atoms (a truncated side chain in a crystal structure, an
#   atom deleted in the Builder) get them rebuilt first -- and the user is
#   told. Before this they were only reported as "not addable".
#   Method (as PDBFixer/Modeller do): the ideal amino acid of the Builder's
#   structure library (structures/13_amino_acids_standard, PDB atom names,
#   3D) is superimposed (Kabsch) on the residue's heavy atoms that exist --
#   at least N, CA, C -- and the missing atoms take its coordinates. When
#   the whole side chain beyond CB is missing, chi1/chi2 are chosen among
#   the staggered values by the fewest contacts with the rest of the system.
#   Bonds (with order/aromaticity) come from the pDynamo component library.
_AMINO_ACID_FILES = {
    'ALA': 'alanine_ala', 'ARG': 'arginine_arg', 'ASN': 'asparagine_asn', 'ASP': 'aspartic_acid_asp',
    'CYS': 'cysteine_cys', 'GLN': 'glutamine_gln', 'GLU': 'glutamic_acid_glu', 'GLY': 'glycine_gly',
    'HIS': 'histidine_his', 'ILE': 'isoleucine_ile', 'LEU': 'leucine_leu', 'LYS': 'lysine_lys',
    'MET': 'methionine_met', 'PHE': 'phenylalanine_phe', 'PRO': 'proline_pro', 'SER': 'serine_ser',
    'THR': 'threonine_thr', 'TRP': 'tryptophan_trp', 'TYR': 'tyrosine_tyr', 'VAL': 'valine_val',
}
_BACKBONE_NAMES = frozenset(('N', 'CA', 'C', 'O', 'OXT'))
_IDEAL_RESIDUES = {}
CLASH_DISTANCE = 3.0          # A, heavy-atom contact penalised when choosing the rotamer
_CHI_VALUES = (-60.0, 180.0, 60.0)
_CHI2_AROMATIC_VALUES = (90.0, -90.0, 0.0, 180.0)


def _ideal_heavy_coordinates(resName):
    """ {atom name: xyz} of the heavy atoms of the ideal amino acid, or None. """
    import numpy as np
    resName = (resName or '').upper()
    if resName in _IDEAL_RESIDUES:
        return _IDEAL_RESIDUES[resName]
    coordinates = None
    filename = _AMINO_ACID_FILES.get(resName)
    if filename is not None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gui', 'windows', 'builder', 'structures',
                            '13_amino_acids_standard', filename + '.mol2')
        try:
            coordinates, section = {}, None
            with open(path) as handle:
                for line in handle:
                    if line.startswith('@<TRIPOS>'):
                        section = line.strip()
                        continue
                    fields = line.split()
                    if section == '@<TRIPOS>ATOM' and len(fields) >= 6 and fields[5].split('.')[0] != 'H':
                        coordinates[fields[1]] = np.array([float(fields[2]), float(fields[3]), float(fields[4])])
        except (OSError, ValueError):
            coordinates = None
    _IDEAL_RESIDUES[resName] = coordinates
    return coordinates


def _template_heavy_bonds(resName):
    """ [(name1, name2, bondType, isAromatic)] between heavy atoms of the component. """
    component = _get_library().GetComponent(resName)
    if component is None:
        return []
    atomic = {a.label: a.atomicNumber for a in (component.atoms or [])}
    return [(b.atomLabel1, b.atomLabel2, b.bondType, bool(b.isAromatic)) for b in (component.bonds or [])
            if atomic.get(b.atomLabel1, 1) != 1 and atomic.get(b.atomLabel2, 1) != 1]


def missing_heavy_atoms(resName, existing_names):
    """ Heavy atoms of the amino-acid template absent from the residue
        (terminus-only OXT excluded); [] for non-amino acids. """
    if (resName or '').upper() not in _AMINO_ACID_FILES:
        return []
    atoms, _bonds, _ = _template_atoms_and_bonds(resName, None)
    if atoms is None:
        return []
    existing = set(existing_names)
    return [name for name, z in atoms.items() if z != 1 and name not in existing]


def _kabsch(P, Q):
    """ (R, t) with Q ~ P @ R + t (rows are points). """
    import numpy as np
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = (Vt.T @ D @ U.T).T
    return R, qc - pc @ R


def _dihedral(p0, p1, p2, p3):
    import numpy as np
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def _rotate_about(points, origin, axis, angle_degrees):
    """ Rodrigues rotation of points (n x 3) about the axis through origin. """
    import numpy as np
    k = axis / np.linalg.norm(axis)
    a = np.radians(angle_degrees)
    p = points - origin
    rotated = p * np.cos(a) + np.cross(k, p) * np.sin(a) + np.outer(p @ k, k) * (1.0 - np.cos(a))
    return rotated + origin


def _beyond(start, blocked, adjacency):
    """ Atoms reachable from `start` without going through `blocked`. """
    seen, stack = {start}, [start]
    while stack:
        for other in adjacency.get(stack.pop(), ()):
            if other not in seen and other != blocked:
                seen.add(other)
                stack.append(other)
    return seen


def rebuild_side_chain(resName, present, environment=None):
    """ {missing heavy atom name: xyz} for one amino acid. present: {name:
        xyz} of its heavy atoms that exist; environment: (n, 3) array of the
        other atoms nearby (rotamer choice). None when it cannot be done
        (no ideal template, or N/CA/C missing). """
    import numpy as np
    ideal = _ideal_heavy_coordinates(resName)
    if ideal is None or not {'N', 'CA', 'C'} <= set(present):
        return None
    missing = [n for n in missing_heavy_atoms(resName, present.keys()) if n in ideal]
    if not missing:
        return {}
    adjacency = {}
    for (a, b, _type, _arom) in _template_heavy_bonds(resName):
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)

    # . internal-coordinate placement: each missing atom X from three known
    #   atoms bonded in a chain (parent P, grandparent G, and D) -- the local
    #   fragment D-G-P of the template is superimposed and X copied, so bond,
    #   angle and dihedral are exactly the template's. (Superimposing the
    #   whole residue at once distorted bonds when the template conformation
    #   differs: ILE CG1-CD1 came out 2.9 A.)
    known = {n: np.asarray(xyz, dtype=float) for n, xyz in present.items() if n in ideal}
    placed = {}
    pending = list(missing)
    while pending:
        progress = False
        for x in list(pending):
            frame = None
            for p in sorted(adjacency.get(x, ())):
                if p not in known:
                    continue
                for g in sorted(adjacency.get(p, ())):
                    if g == x or g not in known:
                        continue
                    d = next((n for n in sorted(adjacency.get(g, ())) if n not in (p, x) and n in known), None)
                    if d is None:
                        d = next((n for n in sorted(adjacency.get(p, ())) if n not in (g, x) and n in known), None)
                    if d is not None:
                        frame = (d, g, p)
                        break
                if frame:
                    break
            if frame is None:
                continue
            R, t = _kabsch(np.array([ideal[n] for n in frame]), np.array([known[n] for n in frame]))
            known[x] = placed[x] = ideal[x] @ R + t
            pending.remove(x)
            progress = True
        if not progress:
            break
    if pending:
        return None

    # . rotamer: only when everything beyond CB is ours to place
    side = [n for n in ideal if n not in _BACKBONE_NAMES and n != 'CB']
    coords = dict(present)
    coords.update(placed)
    if side and 'CB' in coords and all(n in placed for n in side) and environment is not None and len(environment):
        env = np.asarray(environment, dtype=float)
        env = env[np.linalg.norm(env - coords['CA'], axis=1) < 16.0]
        gammas = [n for n in adjacency.get('CB', ()) if n not in ('CA',)]
        gamma = next((g for g in sorted(gammas) if len(adjacency.get(g, ())) > 1), sorted(gammas)[0] if gammas else None)
        if gamma is not None:
            beyond_cb = [n for n in _beyond(gamma, 'CB', adjacency) if n in placed] + \
                        [g for g in gammas if g != gamma and g in placed]
            beyond_cb = list(dict.fromkeys(beyond_cb))
            deltas = [n for n in adjacency.get(gamma, ()) if n != 'CB']
            delta = sorted(deltas)[0] if deltas else None
            beyond_g = [n for n in _beyond(delta, gamma, adjacency) if n in placed] if delta else []
            aromatic = resName.upper() in ('PHE', 'TYR', 'TRP', 'HIS')
            chi2_values = (_CHI2_AROMATIC_VALUES if aromatic else _CHI_VALUES) if delta else (None,)

            def score(trial):
                xyz = np.array([trial[n] for n in trial if n in placed])
                if not len(env):
                    return 0.0
                d = np.linalg.norm(xyz[:, None, :] - env[None, :, :], axis=2)
                return float(np.sum(np.clip(CLASH_DISTANCE - d, 0.0, None) ** 2))

            best = None
            for chi1 in _CHI_VALUES:
                trial = dict(coords)
                current = _dihedral(trial['N'], trial['CA'], trial['CB'], trial[gamma])
                moved = _rotate_about(np.array([trial[n] for n in beyond_cb]), trial['CA'],
                                      trial['CB'] - trial['CA'], chi1 - current)
                for n, xyz in zip(beyond_cb, moved): trial[n] = xyz
                for chi2 in chi2_values:
                    trial2 = dict(trial)
                    if chi2 is not None and beyond_g:
                        current2 = _dihedral(trial2['CA'], trial2['CB'], trial2[gamma], trial2[delta])
                        moved2 = _rotate_about(np.array([trial2[n] for n in beyond_g]), trial2['CB'],
                                               trial2[gamma] - trial2['CB'], chi2 - current2)
                        for n, xyz in zip(beyond_g, moved2): trial2[n] = xyz
                    s = score(trial2)
                    if best is None or s < best[0] - 1e-9:
                        best = (s, trial2)
            coords = best[1]
    return {n: np.asarray(coords[n], dtype=float) for n in missing}


def add_missing_heavy_atoms(old_system, residue_keys):
    """ A new System with the missing heavy atoms of the amino acids in
        `residue_keys` rebuilt (rebuild_side_chain()), same Atom objects and
        coordinates for everything else. Returns ( new system or None,
        { key: [rebuilt names] }, { key: [names that could not be rebuilt] } ). """
    import numpy as np
    from pMolecule import Atom, System
    from pMolecule.Sequence import Sequence
    from pScientific.Geometry3 import Coordinates3

    _ensure_connectivity(old_system)
    residues, _ = build_residue_index(old_system)
    residue_order = [key for key, _names in residues]
    residues_map = build_residue_atom_objects(old_system)
    old_coords = old_system.coordinates3
    all_xyz = np.array([[old_coords[a.index, 0], old_coords[a.index, 1], old_coords[a.index, 2]] for a in old_system.atoms])

    additions, rebuilt, failed = {}, {}, {}
    for key in residue_keys:
        atoms_in_residue = residues_map.get(key) or []
        chain, resName, resSeq = key
        names = [a.label for a in atoms_in_residue]
        missing = missing_heavy_atoms(resName, names)
        if not missing:
            continue
        present = {a.label: all_xyz[a.index] for a in atoms_in_residue if a.atomicNumber != 1}
        mine = set(a.index for a in atoms_in_residue)
        environment = all_xyz[[i for i in range(len(all_xyz)) if i not in mine]] if len(all_xyz) > len(mine) else None
        placed = rebuild_side_chain(resName, present, environment)
        if not placed:
            failed[key] = list(missing)
            continue
        atomic = {a.label: a.atomicNumber for a in (_get_library().GetComponent(resName).atoms or [])}
        additions[key] = [(name, atomic.get(name, 6), placed[name]) for name in missing if name in placed]
        rebuilt[key] = [a[0] for a in additions[key]]
        left = [n for n in missing if n not in placed]
        if left:
            failed[key] = left
    if not additions:
        return None, rebuilt, failed

    captured = {id(a): (old_coords[a.index, 0], old_coords[a.index, 1], old_coords[a.index, 2]) for a in old_system.atoms}
    atoms_list, atom_paths, xyz_list, new_bonds = [], [], [], []
    for key in residue_order:
        chain, resName, resSeq = key
        by_name = {}
        for atom in residues_map.get(key, []):
            atoms_list.append(atom)
            atom_paths.append('{}:{}.{}:{}'.format(chain, resName, resSeq, atom.label))
            xyz_list.append(captured[id(atom)])
            by_name[atom.label] = atom
        new_names = set()
        for (name, atomic_number, xyz) in additions.get(key, []):
            new_atom = Atom.WithOptions(atomicNumber=atomic_number, label=name)
            atoms_list.append(new_atom)
            atom_paths.append('{}:{}.{}:{}'.format(chain, resName, resSeq, name))
            xyz_list.append(tuple(float(v) for v in xyz))
            by_name[name] = new_atom
            new_names.add(name)
        if new_names:
            for (a, b, bond_type, aromatic) in _template_heavy_bonds(resName):
                if (a in new_names or b in new_names) and a in by_name and b in by_name:
                    new_bonds.append((by_name[a], by_name[b], bond_type, aromatic))

    sequence = Sequence.FromAtomPaths(atom_paths, atoms=atoms_list)
    position = {id(atom): k for k, atom in enumerate(sequence.GatherAtoms())}
    all_bonds = [(position[id(b.node1)], position[id(b.node2)], b.type, b.isAromatic) for b in old_system.connectivity.bonds]
    all_bonds += [(position[id(a)], position[id(b)], bond_type, aromatic) for (a, b, bond_type, aromatic) in new_bonds]
    new_system = System.FromSequence(sequence, bonds=all_bonds)
    new_system.label = getattr(old_system, 'label', None)
    coords = Coordinates3.WithExtent(len(atoms_list))
    for i, atom in enumerate(atoms_list):
        coords[atom.index, 0], coords[atom.index, 1], coords[atom.index, 2] = xyz_list[i]
    new_system.coordinates3 = coords
    return new_system, rebuilt, failed


def rebuild_system_with_added_hydrogens(old_system, residue_variant_choices, random_seed=None, rebuild_side_chains=True):
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

    # . [EN] 2026-10-18: amino acids with missing heavy atoms get them rebuilt
    #   first (add_missing_heavy_atoms()), then the hydrogens of the completed
    #   residues are added by the normal path below.
    if rebuild_side_chains:
        standard_keys = [key for key, choice in residue_variant_choices.items()
                         if not (isinstance(choice, tuple) and choice and choice[0] in (ORGANIC_CHOICE_KIND, ORGANIC_PLAN_KIND))]
        heavy_system, rebuilt, not_rebuilt = add_missing_heavy_atoms(old_system, standard_keys)
        if heavy_system is not None:
            report = rebuild_system_with_added_hydrogens(heavy_system, residue_variant_choices, random_seed=random_seed,
                                                         rebuild_side_chains=False)
            if report['system'] is None:
                report['system'] = heavy_system
            for key, names in rebuilt.items():
                report['added'][key] = list(names) + list(report['added'].get(key, []))
            report['rebuilt_side_chains'] = rebuilt
            report['not_rebuilt_side_chains'] = not_rebuilt
            return report

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

        if isinstance(variant_choice, tuple) and variant_choice[0] == ORGANIC_PLAN_KIND:
            plan = variant_choice[1]
            changes = []
            for atom, charge in plan['charges']:
                report.setdefault('old_formal_charges', []).append((atom, getattr(atom, 'formalCharge', 0)))
                atom.formalCharge = charge
                changes.append((atom.label, charge))
            if changes:
                report['charge_changes'][key] = changes
            atoms_to_remove.update(id(a) for a in plan['remove'])
            if plan['remove']:
                report.setdefault('removed', {})[key] = [a.label for a in plan['remove']]
            additions = _name_new_hydrogens(atoms_in_residue, plan['add'], removed=plan['remove'])
            if additions:
                residue_additions[key] = additions
                report['added'][key] = [a[0] for a in additions]
            continue

        if isinstance(variant_choice, tuple) and variant_choice[0] == ORGANIC_CHOICE_KIND:
            ph = variant_choice[1]
            report.setdefault('old_formal_charges', []).extend(
                (a, getattr(a, 'formalCharge', 0)) for a in atoms_in_residue)
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
        heavy_counts = _existing_heavy_bond_counts_for_residue(old_system, atoms_in_residue)
        analysis = analyze_residue(resName, existing_names, variant_label=variant_choice, existing_h_bond_counts=h_counts,
                                   existing_heavy_bond_counts=heavy_counts)
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

    # . [EN] 2026-10-18: bond orders (and formal charges) of every ligand re-
    #   perceived from its final hydrogens -- the PDB bonds are all single and
    #   a proton moved in the review changes the Kekule structure.
    ligand_orders = {}
    for key, variant_choice in residue_variant_choices.items():
        if not (isinstance(variant_choice, tuple) and variant_choice[0] in (ORGANIC_CHOICE_KIND, ORGANIC_PLAN_KIND)):
            continue
        if key not in residue_additions and not any(id(a) in atoms_to_remove for a in residues_map.get(key, [])):
            continue
        atoms_in_residue = residues_map.get(key, [])
        heavy = [a for a in atoms_in_residue if a.atomicNumber != 1]
        adjacency = _residue_adjacency(old_system, atoms_in_residue)
        h_parents = [n for a in atoms_in_residue if a.atomicNumber == 1 and id(a) not in atoms_to_remove
                     for n in adjacency[id(a)][:1]]
        h_parents += [parent for (_name, _z, parent, _xyz) in residue_additions.get(key, [])]
        perceived = perceive_ligand_bond_orders(heavy, h_parents, _residue_heavy_pairs(old_system, heavy))
        if perceived is None:
            continue
        orders, charges = perceived
        ligand_orders.update(orders)
        changes = []
        for atom in heavy:
            charge = charges[id(atom)]
            if charge != int(getattr(atom, 'formalCharge', 0) or 0):
                report.setdefault('old_formal_charges', []).append((atom, getattr(atom, 'formalCharge', 0)))
                atom.formalCharge = charge
                changes.append((atom.label, charge))
        if changes:
            merged = dict(report['charge_changes'].get(key, []))
            merged.update(changes)
            report['charge_changes'][key] = list(merged.items())

    if not residue_additions and not atoms_to_remove:
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
    # [EN] 2026-10-02 (performance): hand pDynamo (i, j, type, isAromatic)
    # tuples with integer positions instead of Bond objects / atom tuples --
    # Bond.FromIterable tests `node in nodes` on a plain list for those
    # (O(atoms) per bond, seconds on a ~6000-atom protein). The positions
    # are those of the gathered sequence atoms, i.e. the new system's order.
    position = {id(atom): k for k, atom in enumerate(sequence.GatherAtoms())}
    all_bonds = []
    bond_type_by_order = {1: BondType.Single, 2: BondType.Double, 3: BondType.Triple}
    for b in old_bonds:
        perceived = ligand_orders.get((id(b.node1), id(b.node2)))
        if perceived is not None:
            all_bonds.append((position[id(b.node1)], position[id(b.node2)],
                              bond_type_by_order.get(perceived[0], BondType.Single), perceived[1]))
            continue
        all_bonds.append((position[id(b.node1)], position[id(b.node2)], b.type, b.isAromatic))
    for (n1, n2, bond_type) in new_bonds:
        all_bonds.append((position[id(n1)], position[id(n2)], bond_type, False))
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
    heavy_atoms = [a for a in atoms_in_residue if a.atomicNumber != 1]
    existing_h_atoms = [a for a in atoms_in_residue if a.atomicNumber == 1]
    result = _ph_protonation(system, heavy_atoms, ph=ph, obabel_bin=obabel_bin)
    if result is None:
        return None
    charges, h_by_parent = result[:2]

    charge_changes = []
    for atom, new_charge in zip(heavy_atoms, charges):
        if new_charge != getattr(atom, 'formalCharge', 0):
            atom.formalCharge = new_charge
            charge_changes.append((atom.label, new_charge))

    entries = [(atom, xyz) for atom in heavy_atoms for xyz in h_by_parent.get(id(atom), [])]
    return {
        'missing_hydrogens': _name_new_hydrogens(atoms_in_residue, entries),
        'charge_changes': charge_changes,
        'atoms_to_remove': existing_h_atoms,
    }


def _ph_protonation(system, heavy_atoms, ph=7.4, obabel_bin=None):
    """ [EN] 2026-10-18: hydrogens/charges of one H-less ligand at a pH. The
        native perception (pdynamo/ligand_perception.py: geometric H counts,
        Wang & Case orders, group pKa rules) by default, as the system
        preparation does since 2026-10-14; OpenBabel only with
        EASYHYBRID_LIGAND_CHEMISTRY=obabel. User's report: OpenBabel 3.1.1 fails
        to kekulize XX01D2VNS2BA-a1291_lig.pdb (a 6-chloro-dihydroquinazoline)
        and makes the benzene ring non-aromatic with an sp3 CH.
        Returns None or ( charges, {id(heavy): [xyz]}, {id(heavy): decision text} ). """
    if os.environ.get('EASYHYBRID_LIGAND_CHEMISTRY', 'native').lower() == 'obabel':
        result = _obabel_protonation(system, heavy_atoms, ph=ph, obabel_bin=obabel_bin)
        return None if result is None else (result[0], result[1], {})
    return _native_protonation(system, heavy_atoms, ph=ph)


def _residue_heavy_pairs(system, heavy_atoms):
    """ [(i, j)] positions in heavy_atoms of the bonds between them. """
    position = {id(a): k for k, a in enumerate(heavy_atoms)}
    _ensure_connectivity(system)
    pairs = set()
    for b in system.connectivity.bonds:
        i, j = position.get(id(b.node1)), position.get(id(b.node2))
        if i is not None and j is not None and i != j:
            pairs.add((min(i, j), max(i, j)))
    return sorted(pairs)


def _native_protonation(system, heavy_atoms, ph=7.4):
    from pScientific import PeriodicTable
    from pdynamo import ligand_perception as LP
    if not heavy_atoms:
        return None
    coords = system.coordinates3
    symbols = [PeriodicTable.Symbol(a.atomicNumber) for a in heavy_atoms]
    xyz = [(coords[a.index, 0], coords[a.index, 1], coords[a.index, 2]) for a in heavy_atoms]
    pairs = _residue_heavy_pairs(system, heavy_atoms)
    try:
        result = LP.perceive_ligand(symbols, xyz, bonds=pairs or None, ph=ph)
    except Exception:
        import traceback
        traceback.print_exc()
        return None
    n = len(heavy_atoms)
    h_by_parent = {}
    for heavy, k in result['added_hydrogens']:
        if _atomic_number_is_hydrogen_free_metal(heavy_atoms[heavy].atomicNumber):
            continue
        h_by_parent.setdefault(id(heavy_atoms[heavy]), []).append(tuple(float(v) for v in result['coords'][k]))
    decisions = {}
    for group, atom, pka, action in result.get('decisions', []):
        if atom < n:
            decisions[id(heavy_atoms[atom])] = '{} ({}{})'.format(
                action, group, '' if pka is None else ', pKa ~{:.1f}'.format(pka))
    return [int(q) for q in result['formal_charges'][:n]], h_by_parent, decisions


def perceive_ligand_bond_orders(heavy_atoms, h_parents, pairs):
    """ Bond orders of a ligand whose hydrogens are final: heavy_atoms, the
        parent (Atom) of each hydrogen, heavy pairs [(i, j)] in heavy_atoms.
        Returns ( {(id a, id b): (order, aromatic)} both orientations,
        {id(heavy): formal charge} ) or None. """
    from pScientific import PeriodicTable
    from pdynamo import ligand_perception as LP
    position = {id(a): k for k, a in enumerate(heavy_atoms)}
    symbols = [PeriodicTable.Symbol(a.atomicNumber) for a in heavy_atoms]
    n = len(symbols)
    all_pairs = list(pairs)
    for parent in h_parents:
        if id(parent) in position:
            symbols.append('H')
            all_pairs.append((position[id(parent)], len(symbols) - 1))
    try:
        result = LP.perceive_ligand(symbols, [(0.0, 0.0, 0.0)] * len(symbols), bonds=all_pairs, protonate_groups=False)
    except Exception:
        import traceback
        traceback.print_exc()
        return None
    aromatic = {(min(i, j), max(i, j)) for (i, j) in result['aromatic']}
    orders = {}
    for i, j, order in result['bonds']:
        if i < n and j < n:
            value = (int(order), (min(i, j), max(i, j)) in aromatic)
            orders[(id(heavy_atoms[i]), id(heavy_atoms[j]))] = value
            orders[(id(heavy_atoms[j]), id(heavy_atoms[i]))] = value
    charges = {id(heavy_atoms[k]): int(result['formal_charges'][k]) for k in range(n)}
    return orders, charges


def _obabel_protonation(system, heavy_atoms, ph=7.4, obabel_bin=None):
    """ OpenBabel's `-p <ph>` on the heavy atoms of one residue, mapped back
        by position. Changes nothing. Returns None (no OpenBabel, failure)
        or ( [formal charge per heavy atom], {id(heavy atom): [(x,y,z), ...]} ). """
    from pScientific import PeriodicTable

    obabel_bin = obabel_bin or shutil.which('obabel')
    if not obabel_bin or not heavy_atoms:
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

    charges = [out_atoms[i]['formal_charge'] for i in range(n_heavy)]
    new_h_by_parent = {}
    for out_atom in out_atoms[n_heavy:]:
        if out_atom['element'] != 'H':
            continue
        parents = [p for p in conect.get(out_atom['serial'], []) if p <= n_heavy]
        if not parents:
            continue
        parent_atom = heavy_atoms[parents[0] - 1]
        if _atomic_number_is_hydrogen_free_metal(parent_atom.atomicNumber):
            continue        # no hydrogens on metals (e.g. Fe of a heme), see is_hydrogen_free_metal()
        new_h_by_parent.setdefault(id(parent_atom), []).append((out_atom['x'], out_atom['y'], out_atom['z']))
    return charges, new_h_by_parent


def _name_new_hydrogens(atoms_in_residue, entries, removed=()):
    """ [(parent_atom, xyz or None)] -> [(name, 1, parent_atom, xyz)], names
        unique in the residue (minus the `removed` atoms). """
    from collections import OrderedDict
    new_h_by_parent = OrderedDict()
    for parent_atom, xyz in entries:
        new_h_by_parent.setdefault(id(parent_atom), []).append((parent_atom, xyz))
    removed_ids = {id(a) for a in removed}

    # . Name new H's -- ALWAYS <= 4 characters (see this function's own
    # docstring for why this is a hard requirement, not cosmetic).
    # Preferred: "H" + parent's own label, when a parent gets exactly
    # one new H AND that name is short enough to fit. Fallback: a plain
    # "H"+running-index, unique within the residue -- used whenever the
    # preferred form would exceed 4 characters (a long parent label, or
    # more than one new H on the same parent, which would otherwise need
    # a "_1"/"_2" suffix that rarely fits).
    existing_names = {a.label for a in atoms_in_residue if id(a) not in removed_ids}
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
        for (parent_atom, xyz) in entries:
            name = None
            if not multiple:
                candidate = 'H' + parent_atom.label
                if len(candidate) <= 4 and candidate not in used_names:
                    name = candidate
            if name is None:
                name = _next_fallback_name()
            used_names.add(name)
            missing_hydrogens.append((name, 1, parent_atom, xyz))
    return missing_hydrogens


# =============================================================================
#  Reviewing / editing the protonation of a ligand -- 2026-10-18
# =============================================================================
# . [EN] User's request: inspect the protonation a ligand gets and change it
#   where needed. A BASELINE is either the ligand's current hydrogens or
#   OpenBabel's prediction at a pH (protonate_non_standard_residue()); the
#   user then adds/removes one proton on any N/O/S atom (formal charge +1/-1
#   each). ligand_protonation_plan() turns that into what
#   rebuild_system_with_added_hydrogens() applies (choice
#   (ORGANIC_PLAN_KIND, plan)). Removed protons are the parent's last H; an
#   added one gets its position from BuildHydrogenCoordinates3FromConnectivity.
ORGANIC_PLAN_KIND = 'organic_plan'
_SITE_ELEMENTS = {7: 'N', 8: 'O', 16: 'S'}
_SITE_MAX_NEIGHBOURS = {7: 4, 8: 2, 16: 2}
_SITE_MAX_CHARGE = {7: 1, 8: 0, 16: 0}       # no oxonium/sulfonium states offered


def _residue_adjacency(system, atoms_in_residue):
    """ {id(atom): [bonded atoms]} for the atoms of one residue (bonds to other
        residues included -- they count for the valence). """
    _ensure_connectivity(system)
    ids = {id(a) for a in atoms_in_residue}
    adjacency = {i: [] for i in ids}
    for b in system.connectivity.bonds:
        a1, a2 = b.node1, b.node2
        if id(a1) in ids:
            adjacency[id(a1)].append(a2)
        if id(a2) in ids:
            adjacency[id(a2)].append(a1)
    return adjacency


def ligand_protonation_baseline(system, atoms_in_residue, ph=None, obabel_bin=None):
    """ ph None: the current hydrogens/formal charges; a number: the pH model
        (_ph_protonation(); None when it fails).
        {'source': 'current'|'ph', 'ph', 'atoms', 'heavy', 'adjacency',
         'h_atoms': {id(heavy): [H atoms now]}, 'h_xyz': {id(heavy): [xyz]} (ph only),
         'n_h': {id(heavy): count}, 'charge': {id(heavy): formal charge}} """
    atoms_in_residue = list(atoms_in_residue)
    heavy = [a for a in atoms_in_residue if a.atomicNumber != 1]
    adjacency = _residue_adjacency(system, atoms_in_residue)
    h_atoms = {id(a): [n for n in adjacency[id(a)] if n.atomicNumber == 1] for a in heavy}
    base = {'source': 'current', 'ph': None, 'decisions': {}, 'atoms': atoms_in_residue, 'heavy': heavy,
            'adjacency': adjacency, 'h_atoms': h_atoms, 'h_xyz': None}
    if ph is None:
        base['n_h'] = {i: len(v) for i, v in h_atoms.items()}
        base['charge'] = {id(a): int(getattr(a, 'formalCharge', 0) or 0) for a in heavy}
        return base
    result = _ph_protonation(system, heavy, ph=ph, obabel_bin=obabel_bin)
    if result is None:
        return None
    charges, h_by_parent, decisions = result
    base.update(source='ph', ph=float(ph), decisions=decisions, h_xyz={id(a): list(h_by_parent.get(id(a), [])) for a in heavy},
                n_h={id(a): len(h_by_parent.get(id(a), [])) for a in heavy},
                charge={id(a): int(c) for a, c in zip(heavy, charges)})
    return base


def _site_group(atom, adjacency):
    """ A short name for the chemical group of a N/O/S atom, from the bond graph. """
    heavy = [n for n in adjacency.get(id(atom), []) if n.atomicNumber != 1]
    z = atom.atomicNumber
    if z == 8:
        if len(heavy) >= 2:
            return 'ether/ester O'
        if not heavy:
            return 'water/hydroxide'
        x = heavy[0]
        terminal_o = [n for n in _neighbours_of(x, adjacency) if n.atomicNumber == 8 and
                      sum(1 for m in _neighbours_of(n, adjacency) if m.atomicNumber != 1) == 1]
        n_o = len(terminal_o)
        if x.atomicNumber == 15:
            return 'phosphate O'
        if x.atomicNumber == 16:
            return 'sulfonate/sulfate O' if n_o >= 2 else 'sulfoxide O'
        if x.atomicNumber == 6:
            return 'carboxyl O' if n_o >= 2 else 'hydroxyl/carbonyl O'
        if x.atomicNumber == 7:
            return 'nitro/N-oxide O'
        return 'O'
    if z == 7:
        if any(_in_ring(atom, n, adjacency) for n in heavy):
            return 'ring N'
        for c in heavy:
            if c.atomicNumber == 6 and sum(1 for n in _neighbours_of(c, adjacency) if n.atomicNumber == 7) >= 2:
                return 'amidine/guanidine N'
        for c in heavy:
            if c.atomicNumber in (6, 16) and any(n.atomicNumber == 8 and sum(1 for m in _neighbours_of(n, adjacency)
                                                                          if m.atomicNumber != 1) == 1
                                                 for n in _neighbours_of(c, adjacency)):
                return 'amide/sulfonamide N'
        return 'amine/amide N' if len(heavy) <= 1 else 'N ({} heavy neighbours)'.format(len(heavy))
    if z == 16:
        return 'thiol/thiolate S' if len(heavy) <= 1 else 'thioether S'
    return ''


def _neighbours_of(atom, adjacency):
    return adjacency.get(id(atom), []) if id(atom) in adjacency else []


def _in_ring(atom, neighbour, adjacency):
    """ True when the bond atom-neighbour is in a ring of the residue (a path
        from neighbour back to atom that does not use that bond). """
    seen, stack = {id(atom)}, [neighbour]
    seen.add(id(neighbour))
    while stack:
        current = stack.pop()
        for n in _neighbours_of(current, adjacency):
            if n is atom and current is not neighbour:
                return True
            if id(n) not in seen and n.atomicNumber != 1 and id(n) in adjacency:
                seen.add(id(n))
                stack.append(n)
    return False


def proton_state_text(symbol, n_h, charge):
    """ 'NH3+', 'O-', 'OH', 'N', 'SH' ... """
    text = symbol + ('H' if n_h == 1 else 'H{}'.format(n_h) if n_h > 1 else '')
    if charge:
        text += ('+' if charge > 0 else '-') * min(abs(charge), 1) + (str(abs(charge)) if abs(charge) > 1 else '')
    return text


def ligand_protonation_sites(baseline):
    """ The N/O/S atoms of the ligand, each with the states offered:
        [{'atom', 'label', 'symbol', 'group', 'bonded', 'n_h', 'charge',
          'states': [(delta, text), ...]}] -- delta = protons added (-1, 0, +1;
        resulting formal charge -1..+1 for N, -1..0 for O/S, neighbour count
        within the element's limit). """
    adjacency = baseline['adjacency']
    sites = []
    for atom in baseline['heavy']:
        symbol = _SITE_ELEMENTS.get(atom.atomicNumber)
        if symbol is None:
            continue
        n_h, charge = baseline['n_h'][id(atom)], baseline['charge'][id(atom)]
        n_heavy = sum(1 for n in adjacency[id(atom)] if n.atomicNumber != 1)
        states = []
        for delta in (-1, 0, 1):
            h, q = n_h + delta, charge + delta
            if delta and (h < 0 or q < -1 or q > _SITE_MAX_CHARGE[atom.atomicNumber]
                          or n_heavy + h > _SITE_MAX_NEIGHBOURS[atom.atomicNumber]):
                continue
            states.append((delta, proton_state_text(symbol, h, q)))
        bonded = ', '.join(n.label for n in adjacency[id(atom)] if n.atomicNumber != 1)
        sites.append({'atom': atom, 'label': atom.label, 'symbol': symbol, 'group': _site_group(atom, adjacency),
                      'rule': baseline.get('decisions', {}).get(id(atom), ''),
                      'bonded': bonded, 'n_h': n_h, 'charge': charge, 'states': states})
    return sites


def ligand_net_charge(baseline, deltas=None):
    deltas = deltas or {}
    return sum(baseline['charge'].values()) + sum(deltas.values())


def ligand_protonation_plan(baseline, deltas=None):
    """ {id(heavy atom): protons added (+1) / removed (-1)} on top of the
        baseline -> {'remove': [H atoms], 'add': [(parent, xyz or None)],
        'charges': [(atom, formal charge)], 'net_charge', 'source', 'ph'}. """
    deltas = {k: v for k, v in (deltas or {}).items() if v}
    remove, add, charges = [], [], []
    for atom in baseline['heavy']:
        i = id(atom)
        delta = deltas.get(i, 0)
        if baseline['source'] == 'current':
            current = baseline['h_atoms'][i]
            if delta < 0:
                remove.extend(current[len(current) + delta:])
        else:
            remove.extend(baseline['h_atoms'][i])
            xyz = baseline['h_xyz'][i]
            add.extend((atom, p) for p in xyz[:len(xyz) + min(delta, 0)])
        add.extend((atom, None) for _ in range(max(delta, 0)))
        charge = baseline['charge'][i] + delta
        if charge != int(getattr(atom, 'formalCharge', 0) or 0):
            charges.append((atom, charge))
    return {'remove': remove, 'add': add, 'charges': charges, 'atoms': baseline['atoms'],
            'net_charge': ligand_net_charge(baseline, deltas), 'source': baseline['source'], 'ph': baseline['ph'],
            'n_changes': len(deltas)}


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
