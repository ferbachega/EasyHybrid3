#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Pure-Python PDBQT structural parser
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      Parses PDBQT (AutoDock/Vina) ATOM/HETATM records directly into
#      plain-Python structures, without ever going through OpenBabel or
#      pDynamo. No GTK/vismol imports here (same convention as
#      vina_runner.py/autodock_gpu_runner.py) so this module stays
#      standalone-testable.
#
#      PDBQT is PDB-fixed-column-compatible through occupancy/bfactor,
#      followed by two PDBQT-specific trailing fields (partial charge,
#      AutoDock atom type) that are NOT fixed-width in real Vina/
#      AutoDockTools output -- they are parsed via a plain .split() on
#      the tail of the line instead, confirmed against real files:
#      "...  0.00  0.00    -0.244 NA" -> tail.split() == ['-0.244', 'NA'].
#


def _element_from_name(name):
    """ Derives a chemical element from a PDB/PDBQT atom NAME, ported
        CASE-SENSITIVELY from vismol/model/atom.py's own Atom._get_symbol()
        (the established, battle-tested logic already used everywhere
        else in this codebase for the exact same problem) -- NOT a
        naive .capitalize()-based guess.

        This case-sensitivity is load-bearing, not cosmetic: PDB/PDBQT
        atom NAMES are conventionally ALL-CAPS (e.g. "CA" for a protein
        alpha-carbon, "CB", "CL1"), while real two-letter element
        SYMBOLS are Title-case (e.g. "Ca" calcium, "Cl" chlorine). An
        earlier version of this function called .capitalize() on the
        name before comparing it against a two-letter-element
        whitelist, which THROWS AWAY that case distinction -- "CA"
        became "Ca" and was misidentified as CALCIUM instead of carbon.
        Confirmed as a real, systemic bug against a real receptor
        PDBQT: every single alpha-carbon (atom name "CA", one per
        residue) in the whole protein got built as calcium, inflating
        each one's covalent radius for VismolObject.
        find_bonded_and_nonbonded_atoms()'s distance-based bond
        perception and producing spurious bonds to many nearby atoms --
        this is what a user reported as "the structure shows wrong
        bonds, like they're out of order."

        NEVER derived from the PDBQT "autodock_type" field, which is a
        docking-specific pseudo-element (e.g. 'A' aromatic carbon,
        'OA'/'NA'/'SA' H-bond acceptor variants, 'HD' donor hydrogen)
        that must not be fed into a chemical periodic-table lookup.
    """
    name = ''.join(ch for ch in name.strip() if not ch.isnumeric())
    if len(name) >= 3:
        name = name[:2]
    if not name:
        return 'X'

    first = name[0]
    second = name[1] if len(name) > 1 else ''

    if first == 'H':
        if second == 'g': return 'Hg'
        if second == 'e': return 'He'
        return 'H'
    if first == 'C':
        if second == 'a': return 'Ca'
        if second in ('l', 'L'): return 'Cl'
        if second == 'd': return 'Cd'
        if second in ('u', 'U'): return 'Cu'
        return 'C'
    if first == 'N':
        if second in ('i', 'I'): return 'Ni'
        if second == 'a': return 'Na'
        if second == 'e': return 'Ne'
        if second == 'b': return 'Nb'
        return 'N'
    if first == 'O':
        if second == 's': return 'Os'
        return 'O'
    if first == 'S':
        if second == 'I': return 'Si'
        if second == 'e': return 'Se'
        return 'S'
    if first == 'P':
        if second == 'd': return 'Pd'
        if second == 'b': return 'Pb'
        if second == 'o': return 'Po'
        return 'P'
    if first == 'Z':
        if second == 'r': return 'Zr'
        return 'Zn'
    if first == 'F':
        if second in ('E', 'e'): return 'Fe'
        return 'F'
    if first == 'M':
        if second in ('n', 'N'): return 'Mn'
        if second == 'o': return 'Mo'
        if second == 'G': return 'Mg'
        return 'X'
    if first == 'B':
        if second == 'r': return 'Br'
        return 'B'
    if first == 'I':
        return 'I'
    return 'X'


def parse_pdbqt_atom_line(line):
    """ Parses one ATOM/HETATM PDBQT line into a plain dict:
        {serial, name, resname, chain, resseq, x, y, z, occupancy,
        bfactor, charge, autodock_type, element}.

        Uses standard PDB fixed columns for everything through
        occupancy/bfactor (name 13-16, resName 18-20, chainID 22,
        resSeq 23-26, x 31-38, y 39-46, z 47-54, occupancy 55-60,
        bfactor 61-66, all 1-indexed/inclusive -- same columns this
        codebase already trusts elsewhere, see export_special_PDB() in
        pdynamo/pDynamo2EasyHybrid/helpers.py and
        vina_runner.read_pdb_coordinates()/read_pdb_elements()), then
        `line[66:].split()` for the two trailing PDBQT-only fields
        (charge, autodock_type) since those are NOT fixed-width.
    """
    name = line[12:16].strip()
    tail = line[66:].split()
    charge = float(tail[0]) if len(tail) >= 1 else 0.0
    autodock_type = tail[1] if len(tail) >= 2 else ''
    chain = line[21].strip()
    return {
        'serial': int(line[6:11]),
        'name': name,
        'resname': line[17:20].strip(),
        'chain': chain if chain else ' ',
        'resseq': int(line[22:26]),
        'x': float(line[30:38]),
        'y': float(line[38:46]),
        'z': float(line[46:54]),
        'occupancy': float(line[54:60]) if line[54:60].strip() else 0.0,
        'bfactor': float(line[60:66]) if line[60:66].strip() else 0.0,
        'charge': charge,
        'autodock_type': autodock_type,
        'element': _element_from_name(name),
    }


def _is_atom_record(line):
    return line.startswith('ATOM') or line.startswith('HETATM')


def parse_pdbqt_topology(pdbqt_path):
    """ Parses the FIRST MODEL only (or the whole file if it has no
        MODEL/ENDMDL blocks at all -- e.g. a receptor file) into an
        ordered list of atom dicts (see parse_pdbqt_atom_line()), each
        additionally carrying 'index' = 0-based order of appearance.

        ATOM/HETATM record order is preserved verbatim -- this order IS
        the atom-numbering contract every VObject built from this
        topology depends on. ROOT/ENDROOT/BRANCH/ENDBRANCH/TORSDOF/
        REMARK/TER lines are skipped for topology purposes.
    """
    atoms = []
    in_first_model = True
    with open(pdbqt_path, 'r') as handle:
        for line in handle:
            if line.startswith('ENDMDL'):
                break
            if not in_first_model:
                continue
            if _is_atom_record(line):
                atom = parse_pdbqt_atom_line(line)
                atom['index'] = len(atoms)
                atoms.append(atom)
    return atoms


def parse_pdbqt_models(pdbqt_path):
    """ Splits a multi-MODEL ligand-pose PDBQT (MODEL/ENDMDL blocks)
        into a list of per-model records:
        [{'model': int, 'atoms': [atom_dict, ...]}, ...]

        A single-pose PDBQT with no MODEL/ENDMDL lines at all returns a
        single-element list with model=1. This is a pure-Python
        replacement for vina_runner.extract_pose_as_pdb()'s obabel-
        based MODEL-splitting -- no subprocess call needed.

        Do not trust `model` as a stable pose/mode or GA-run number for
        AutoDock-GPU results -- dlg_to_docked_pdbqt() already renumbers
        MODEL == GA run number by its own convention; callers should
        keep passing mode/run explicitly from their own results dict
        rather than relying on this field alone.
    """
    models = []
    current_model = None
    current_atoms = None
    saw_any_model_line = False
    with open(pdbqt_path, 'r') as handle:
        for line in handle:
            if line.startswith('MODEL'):
                saw_any_model_line = True
                current_model = int(line.split()[1])
                current_atoms = []
            elif line.startswith('ENDMDL'):
                if current_atoms is not None:
                    for i, atom in enumerate(current_atoms):
                        atom['index'] = i
                    models.append({'model': current_model, 'atoms': current_atoms})
                current_model = None
                current_atoms = None
            elif _is_atom_record(line):
                if current_atoms is None:
                    current_atoms = []
                current_atoms.append(parse_pdbqt_atom_line(line))

    if not saw_any_model_line:
        atoms = parse_pdbqt_topology(pdbqt_path)
        return [{'model': 1, 'atoms': atoms}]
    return models


def pdbqt_topology_names_match(topology_a, topology_b):
    """ True iff (name, resname, chain, resseq) matches at every index
        for two topologies of equal length -- the PDBQT-native
        replacement for vina_runner.read_pdb_elements()'s element-only
        check, and strictly stronger (name-level, not just element-
        level) since PDBQT guarantees identical atom order/naming
        across poses of the SAME docking job by construction (unlike
        the OpenBabel-PDB-extraction path this replaces, which did not
        guarantee this -- that mismatch is why the original, weaker
        element-only check existed at all).
    """
    if len(topology_a) != len(topology_b):
        return False
    for a, b in zip(topology_a, topology_b):
        if (a['name'], a['resname'], a['chain'], a['resseq']) != \
           (b['name'], b['resname'], b['chain'], b['resseq']):
            return False
    return True


def drop_duplicate_atom_names(topology):
    """ Drops every atom past the FIRST occurrence of a colliding
        (chain, resname, resseq, name) key within `topology`, keeping
        the surviving atoms in their original relative order and
        re-deriving their 'index' to stay contiguous (0..N-1).

        RECEPTOR policy -- use this for protein/receptor topologies,
        NOT for ligands (see rename_duplicate_atom_names() below for
        why the two need opposite policies).

        A real protein residue never legitimately repeats an atom name
        within itself -- a collision here means the source PDBQT
        concatenated overlapping structural data for that residue, not
        that two genuinely different atoms happen to share a generic
        label. Confirmed against a real receptor file (7CO0's
        flexible_receptor_PDBQT/receptor.pdbqt): MET77 and GLU215 each
        appear TWICE, sharing a near-coincident CA anchor (~0.002-0.13
        A apart) but with the REST of the sidechain diverging
        increasingly further out (up to ~5 A for GLU215's OE1/OE2) --
        i.e. two different, overlapping conformations of the same
        residue's sidechain, most likely a leftover from concatenating
        AutoDockTools' separately-generated rigid-receptor and
        flexible-residue PDBQT outputs into one file without excluding
        the flexible residue from the rigid one.

        Feeding two overlapping/divergent conformations of one residue
        into distance-based bond perception
        (VismolObject.find_bonded_and_nonbonded_atoms(), used for every
        PDBQT-built object -- see vismol_pdbqt_builder.py) creates
        spurious bonds between the two near-coincident copies, which is
        exactly what a user reported as "the structure shows wrong
        bonds, like they're out of order." Dropping the redundant copy
        instead guarantees exactly one, internally self-consistent
        conformation per residue, so there is nothing left for bond
        perception to get confused by. A residue with no collisions is
        left byte-identical.

        Do NOT use this for ligands: a small-molecule PDBQT can
        legitimately have SEVERAL genuinely distinct atoms sharing one
        generic name (confirmed against a real fixture,
        adtVinaFR/7CO0/ligands_PDBQT/orca.finalensemble_002.pdbqt: 3
        real, chemically distinct carbons all named plain "C") --
        dropping there would silently delete real atoms (a 5-atom
        ligand collapsing to 3), which is exactly what a user reported
        next as "the ligand structures aren't opening correctly" right
        after this function's policy was first applied universally.
        rename_duplicate_atom_names() is the correct policy there.
    """
    seen = set()
    deduped = []
    for atom in topology:
        key = (atom['chain'], atom['resname'], atom['resseq'], atom['name'])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(atom)
    topology[:] = deduped
    for i, atom in enumerate(topology):
        atom['index'] = i


def rename_duplicate_atom_names(topology):
    """ Renames atoms that collide on (chain, resname, resseq, name)
        within `topology`, in place, so every atom gets a unique
        (chain, resname, resseq, name) identity, WITHOUT dropping any
        atom -- same strategy as util/pdb_tools.py's
        dedupe_pdb_atom_names(), applied directly to an in-memory
        topology list.

        LIGAND policy -- use this for small-molecule/ligand topologies,
        NOT for receptors (see drop_duplicate_atom_names() above for
        the receptor case, which needs the opposite policy: a real
        protein residue's name collision means redundant/overlapping
        data to be dropped, not genuinely distinct atoms to preserve).

        Confirmed against a real ligand fixture (5-atom
        R-epichlorohydrin PDBQT): 3 of its 5 atoms are all named plain
        "C" -- three REAL, chemically distinct carbons that a
        generic/simple ligand-prep step never gave individual labels
        to. Every atom of a residue that has ANY name collision gets
        renamed to element + a running PER-ELEMENT count across the
        whole residue (e.g. all three carbons become C1/C2/C3, even the
        one that happened to be unique already) -- keying the running
        count per ELEMENT (not per original name) avoids two different
        original names independently producing the same new name (e.g.
        "CA"->"C1"/"C2" and "CB"->"C1"/"C2" colliding with each other).
        This only disambiguates the PATH used to build the pDynamo
        Sequence (Sequence.FromAtomPaths raises on duplicate paths); it
        does not change coordinates, the AutoDock type, or the charge
        on any atom. A residue with no collisions is left byte-identical.
    """
    groups = {}
    for atom in topology:
        key = (atom['chain'], atom['resname'], atom['resseq'])
        groups.setdefault(key, []).append(atom)

    for atoms_in_residue in groups.values():
        name_counts = {}
        for atom in atoms_in_residue:
            name_counts[atom['name']] = name_counts.get(atom['name'], 0) + 1
        if not any(count > 1 for count in name_counts.values()):
            continue
        running = {}
        for atom in atoms_in_residue:
            element = atom['element']
            running[element] = running.get(element, 0) + 1
            atom['name'] = '{}{}'.format(element, running[element])


def get_coords_array(model_record):
    """ model_record: one element of parse_pdbqt_models()'s return, or
        any dict with an 'atoms' key holding a topology list.
        Returns np.ndarray, shape (n_atoms, 3), float32, atom-index order.
    """
    import numpy as np
    atoms = model_record['atoms'] if isinstance(model_record, dict) and 'atoms' in model_record else model_record
    return np.array([[a['x'], a['y'], a['z']] for a in atoms], dtype=np.float32)
