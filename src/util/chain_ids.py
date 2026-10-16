#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  chain_ids.py
#
#  [EN] 2026-09-29 BUG FIX (user report: the cartoon of
#  examples/pkl/citrate_synthase_CHARMM_AM1.pkl -- a truncated structure --
#  connected stretches that should be broken).
#
#  EasyHybrid used the FIRST CHARACTER of each pDynamo sequence entity label
#  as the chain ID (`entityLabel[0:1]`). That is right for PDB-derived
#  systems, whose entities are already called "A", "B", ... -- but CHARMM/
#  PSF segment IDs are 4 characters and commonly share their first letter:
#  citrate synthase has segments AAAA (protein A), AAAW (waters A), AABA
#  (protein B, truncated) and AABW -- all four became chain "A". Residues of
#  the two proteins with the same number were then merged into ONE vismol
#  Residue (ILE 36 got 36 atoms, with the N/CA of one residue and the C of
#  another 11 A away), so every per-residue feature saw garbage: the cartoon
#  bridged the truncation gaps, chain selection grabbed both proteins, etc.
#
#  Rule: the chain ID is the FULL entity label. PDB-derived entities are
#  already single characters ("A", "B", ...), so nothing changes there;
#  multi-character labels (CHARMM/PSF segments "AAAA", "PROA", "WAT1",
#  ligands "OAA1"...) are now kept whole. Entity labels are unique within a
#  pDynamo system, and the full label is exactly what pDynamo's own atom
#  patterns expect ("AAAA:ILE.36:CA" -- the pDynamo Selection window builds
#  its pattern from the vismol chain name; with "A" it matched 0 atoms).
# ============================================================================


def entity_chain_id_map ( entity_labels ):
    """ {entity_label: chain_id} for an iterable of entity labels (see the
        module docstring for the rule). Duplicated labels are fine. Kept as
        a mapping (not inlined) so the rule lives in one place for every
        caller (session.py, hydrogen_builder.py). """
    return { label: label for label in dict.fromkeys ( entity_labels ) }


def system_chain_id_map ( system ):
    """ entity_chain_id_map() over every entity that owns an atom of a
        pDynamo `system` (atom.parent.parent is the atom's entity). """
    labels = [ ]
    for atom in system.atoms:
        try:
            labels.append ( atom.parent.parent.label )
        except AttributeError:
            pass
    return entity_chain_id_map ( labels )


def chain_id_for_atom ( atom, mapping = None ):
    """ Chain ID of a pDynamo atom: mapping[entity label] when given (and
        known), else the historical first-character rule. """
    label = atom.parent.parent.label
    if mapping is not None and label in mapping:
        return mapping[label]
    return label[0:1]
