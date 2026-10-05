#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  builder_solvation.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  [EN] 2026-10-03: applies solvation.solvate() to the Builder's target
#  object. Force-field independent: waters become residues HOH (atoms O, H1,
#  H2, two O-H bonds) in their own chain, ions become one-atom residues (NA,
#  K, CL, ... with their formal charge) in another chain -- the names every
#  force field / PDB tool expects. [EN] 2026-10-03: any solvent of the
#  library (solvent_library.py): its residue name, atom names and bond orders
#  come from its .mol2 file; chain "W" for water, "S" for other solvents. Done in ONE rebuild through the same
#  machinery as Undo (atom_ops._snapshot_builder_state() +
#  _restore_builder_state()), so a protein + thousands of waters costs one
#  batched rebuild instead of one add_atom() redraw per atom, and a single
#  Undo removes the whole solvation (cell included).
# ============================================================================
from gui.windows.builder import solvation
from gui.windows.builder.solvation import SolvationError


def _unique_chain_name ( vismol_object, preferred ):
    taken = set ( str ( name ) for name in ( getattr ( vismol_object, "chains", None ) or { } ) )
    if preferred not in taken:
        return preferred
    k = 2
    while "{}{}".format ( preferred, k ) in taken:
        k += 1
    return "{}{}".format ( preferred, k )


def solute_charge ( vismol_object ):
    """ Total formal charge of the object (the Builder's own formal charges:
        protonation tools, residue templates, fragments). """
    return int ( sum ( int ( getattr ( atom, "formal_charge", 0 ) or 0 )
                       for atom in vismol_object.atoms.values ( ) ) )


def solvate_builder_object ( vismol_object, water_chain = "W", ion_chain = "I", **parameters ):
    """ Solvates `vismol_object` in place (see solvation.solvate() for the
        parameters). The solute is moved so the box starts at the origin
        (that is where vismol draws cells). Returns the SolvationResult.
        Raises SolvationError with a user-facing message. """
    from gui.windows.builder.atom_ops import ( push_undo_snapshot, _snapshot_builder_state,
                                               _restore_builder_state )
    if getattr ( vismol_object, "builder_cell", None ) or getattr ( vismol_object, "builder_solvation", None ):
        raise SolvationError ( "This molecule was already solvated in the Builder. "
                               "Undo the previous solvation first." )

    snapshot = _snapshot_builder_state ( vismol_object )
    atoms    = snapshot["atoms"]
    import numpy as np
    solute_xyz = np.array ( [ [ a["x"], a["y"], a["z"] ] for a in atoms ], dtype = np.float64 ).reshape ( -1, 3 )
    charge     = solute_charge ( vismol_object )

    result = solvation.solvate ( solute_xyz, solute_charge = charge, **parameters )

    # solute moved into the middle of the box (and rotated, if aligned);
    # a sphere leaves it where it is
    for a, xyz in zip ( atoms, result.solute_xyz ):
        a["x"], a["y"], a["z"] = float ( xyz[0] ), float ( xyz[1] ), float ( xyz[2] )

    # solvent molecules: atoms/names/residue/bond orders from the solvent's
    # library file (solvent_library.py); chain "W" for water, "S" otherwise
    solvent = result.solvent
    solvent_chain = _unique_chain_name ( vismol_object, water_chain if solvent.is_water else "S" )
    ion_chain     = _unique_chain_name ( vismol_object, ion_chain )
    bonds  = set ( snapshot["manual_bonds"] )
    orders = dict ( snapshot["manual_bond_orders"] )
    n_at   = solvent.n_atoms
    first  = len ( atoms )
    for k, molecule in enumerate ( result.molecules ):
        for name, symbol, xyz in zip ( solvent.names, solvent.elements, molecule ):
            atoms.append ( { "symbol": symbol, "name": name, "chain_id": solvent_chain, "resi": k + 1,
                             "resn": solvent.residue, "x": float ( xyz[0] ), "y": float ( xyz[1] ),
                             "z": float ( xyz[2] ), "formal_charge": 0 } )
        base = first + k * n_at
        for ( i, j, order ) in solvent.bonds:
            pair = ( base + min ( i, j ), base + max ( i, j ) )
            bonds.add ( pair )
            if order != 1:
                orders[pair] = order
    for k, ( label, element, name, charge_k, xyz ) in enumerate ( result.ions ):
        atoms.append ( { "symbol": element, "name": name, "chain_id": ion_chain, "resi": k + 1,
                         "resn": name, "x": float ( xyz[0] ), "y": float ( xyz[1] ), "z": float ( xyz[2] ),
                         "formal_charge": int ( charge_k ) } )
    snapshot["manual_bonds"]       = bonds
    snapshot["manual_bond_orders"] = orders
    snapshot["builder_cell"]       = result.cell
    info = { "shape": result.region.shape, "solvent": solvent.name, "molecules": result.n_waters,
             "ions": len ( result.ions ) }
    if result.sphere is not None:
        info["centre"] = [ float ( v ) for v in result.sphere[0] ]
        info["radius"] = float ( result.sphere[1] )
    snapshot["builder_solvation"]  = info

    push_undo_snapshot ( vismol_object )
    _restore_builder_state ( vismol_object, snapshot )
    return result
