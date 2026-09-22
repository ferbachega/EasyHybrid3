#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  atom_types.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Live, on-demand DYFF atom typing for a selection of atoms, plus a
#  per-atom manual override mechanism -- the user's own explicit request:
#  "para uma selecao de atomos dentro do modo de edicao, o usuario
#  pudesse verificar qual tipo de atomos foi associado, e mudar isso, se
#  for necessario ... quero verificar o tipo ainda no builder".
#
#  Deliberately does NOT require the system to already have DYFF
#  assigned (p_session.define_MMModel(), see project_dyff_force_field_
#  assignment) -- runs the SAME 3-call sequence MMModelDYFF.BuildModel()
#  itself uses (DetermineAtomGeometry -> DetermineAtomOxidationState ->
#  MMAtomTyper.TypeAtoms(), confirmed by reading pMolecule/MMModel/
#  MMModel.py directly -- skipping the first two gives WRONG types, this
#  project's OWN earlier "sp2 carbon looks sp3" false-alarm this session
#  was caused by exactly that mistake in a standalone test script) on a
#  throwaway pDynamo System built fresh from the CURRENT Builder object
#  state (empty_object._build_pdynamo_system_from_vismol_object(), same
#  helper the Builder's own sync already uses) -- so the Atom Types
#  window works at ANY point while editing, not just after a real DYFF
#  assignment.
# ============================================================================
import os


def _dyff_type_labels_by_atomic_number ( parameter_set = "dyff-1.0" ):
    """ [EN] {atomic_number: sorted [label, ...]} for every DYFF atom
    type this install's parameter set actually defines -- read directly
    from the SAME MMAtomTyper/MMAtomTypeContainer machinery DYFF's own
    typer uses (pMolecule.MMModel.MMAtomTyper.FromPath(...).mmAtomTypes.
    rawItems, a dict keyed by label -> MMAtomType with its own
    .atomicNumber/.label) -- not a hand-copied duplicate of atomTypes.
    yaml's contents, so this can never drift out of sync with whatever
    parameter set is actually installed. Used both for the Atom Types
    window's per-row tooltip ("valid types for this atom's element") and
    to VALIDATE a manually-typed override before accepting it.

    Cached at module level (this is static, load-once data -- the
    parameter set doesn't change while the app is running). """
    cached = _dyff_type_labels_by_atomic_number._cache.get ( parameter_set )
    if cached is not None:
        return cached

    from pMolecule.MMModel.MMAtomTyper import MMAtomTyper
    base_dir = os.path.join ( os.environ["PDYNAMO3_PARAMETERS"], "forceFields", "dyff", parameter_set )
    typer = MMAtomTyper.FromPath ( base_dir )

    by_number = { }
    if typer.mmAtomTypes is not None:
        for label, mm_atom_type in typer.mmAtomTypes.rawItems.items ( ):
            by_number.setdefault ( mm_atom_type.atomicNumber, [ ] ).append ( label )
    for labels in by_number.values ( ):
        labels.sort ( )

    _dyff_type_labels_by_atomic_number._cache[parameter_set] = by_number
    return by_number

_dyff_type_labels_by_atomic_number._cache = { }


def valid_type_labels_for_symbol ( vm_session, symbol, parameter_set = "dyff-1.0" ):
    """ Sorted list of every valid DYFF type label for this element (by
    symbol, e.g. "C" -> ["C:Lin", "C:Res", "C:Tet", "C:Tri"]) -- empty
    list if this parameter set defines none for that element (or the
    symbol isn't recognised at all). """
    entry = vm_session.periodic_table.elements_by_symbol.get ( symbol )
    if entry is None:
        return [ ]
    atomic_number = entry[0]
    return list ( _dyff_type_labels_by_atomic_number ( parameter_set ).get ( atomic_number, [ ] ) )


def compute_atom_types ( vismol_object, atom_ids = None, parameter_set = "dyff-1.0" ):
    """ Returns { atom_id: (perceived_type, effective_type) } for every
    atom_id in `atom_ids` (or every atom in vismol_object if None) --
    perceived_type is whatever MMAtomTyper.TypeAtoms() itself came up
    with, PURELY from the current element/connectivity/bond-order/
    geometry; effective_type is perceived_type UNLESS vismol_object.
    manual_atom_type_overrides has an entry for that atom_id, in which
    case the override wins -- callers (the Atom Types window) show BOTH
    so it's visually obvious which rows are auto-perceived vs manually
    corrected.

    Returns (None, override_or_None) for an atom_id MMAtomTyper left
    untyped (its own CheckUntypedAtoms() raising a whole-system
    MMModelError is caught HERE and turned into a per-atom None instead
    -- one untyped atom elsewhere in the molecule shouldn't stop the
    user from inspecting/fixing the ones they actually selected). """
    from gui.windows.builder.empty_object import _build_pdynamo_system_from_vismol_object
    from pMolecule.ConnectivityUtilities import DetermineAtomGeometry, DetermineAtomOxidationState
    from pMolecule.MMModel.MMAtomTyper import MMAtomTyper
    from pMolecule.MMModel.MMModelError import MMModelError

    overrides = getattr ( vismol_object, "manual_atom_type_overrides", None ) or { }
    scope = list ( vismol_object.atoms.keys ( ) ) if atom_ids is None else list ( atom_ids )

    system = _build_pdynamo_system_from_vismol_object ( vismol_object, label = vismol_object.name )
    DetermineAtomGeometry ( system.connectivity )
    DetermineAtomOxidationState ( system.connectivity )

    base_dir = os.path.join ( os.environ["PDYNAMO3_PARAMETERS"], "forceFields", "dyff", parameter_set )
    typer = MMAtomTyper.FromPath ( base_dir )
    try:
        ( perceived_types, _atom_charges ) = typer.TypeAtoms ( system.connectivity, system.sequence, None )
    except MMModelError:
        # [EN] Whole-system typing failed (at least one atom untyped
        # anywhere in the molecule, not necessarily among the ones the
        # user actually selected) -- fall back to a per-atom None for
        # EVERY atom rather than surfacing a raw exception here; the
        # window shows "?" for those and the user can still set a
        # manual override on them directly.
        perceived_types = [ None ] * len ( system.connectivity.nodes )

    result = { }
    for atom_id in scope:
        perceived = perceived_types[atom_id] if atom_id < len ( perceived_types ) else None
        override  = overrides.get ( atom_id )
        effective = override if override is not None else perceived
        result[atom_id] = ( perceived, effective )
    return result


def set_manual_atom_type_override ( vismol_object, atom_id, type_label ):
    """ Records (or clears, if type_label is None/empty) a manual DYFF
    type override for a single atom -- same storage convention as
    manual_bonds/manual_bond_orders (atom_ops.py's add_atom() docstring):
    a plain dict, lazily created, keyed by atom_id. Does NOT validate
    type_label against valid_type_labels_for_symbol() itself -- the
    Atom Types window's own cell-edited handler does that BEFORE calling
    this, so a rejected edit never reaches here at all. """
    if not hasattr ( vismol_object, "manual_atom_type_overrides" ) or vismol_object.manual_atom_type_overrides is None:
        vismol_object.manual_atom_type_overrides = { }
    if type_label:
        vismol_object.manual_atom_type_overrides[atom_id] = type_label
    else:
        vismol_object.manual_atom_type_overrides.pop ( atom_id, None )
