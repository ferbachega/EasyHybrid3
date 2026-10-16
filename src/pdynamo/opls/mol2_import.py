#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: fallback .mol2 import
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      pDynamo's MOL2 reader rejects some charged aromatic molecules
#      ("Error converting input connectivity (undefined bonds or implicit
#      hydrogens)": thiazolium, pyridinium, cationic purines, some
#      carbamates) because a mol2 file carries no formal charges and its
#      Kekule assignment then fails. This fallback builds the system the
#      way the Builder does (Connectivity + System.FromConnectivity, no
#      residues): bond orders from the file with aromatic bonds Kekulized
#      by the project's perceiver, valence deficits repaired (C=N+ in
#      cationic rings) and formal charges from the valences -- the same
#      chemistry as the OPLS ligand parametrization (ligand.py).
#

import os

import numpy as np


def system_from_mol2 ( path, label = None ):
    """ pDynamo System (no sequence) from a .mol2 file. """
    from pMolecule import Atom, Bond, BondType, Connectivity, ConvertInputConnectivity, System
    from pScientific.Geometry3 import Coordinates3
    from pdynamo.opls import ligand as L
    atoms, bonds = L._parse_mol2 ( path )
    if not atoms:
        raise ValueError ( "{}: no atoms.".format ( path ) )
    elements = [ a[1] for a in atoms ]
    kekule, aromatic = L._kekulize ( elements, bonds )
    kekule = L._repair_valence_deficits ( elements, kekule, [ a[2] for a in atoms ] )
    charges = L._formal_charges ( elements, kekule )
    # . unique atom names (mol2 files often name every atom by its element)
    counts, names, used = { }, [ ], set ( )
    for ( name, element, _, _ ) in atoms: counts[name] = counts.get ( name, 0 ) + 1
    for ( name, element, _, _ ) in atoms:
        new = name[:4] if counts[name] == 1 else None
        if new is None or new in used:
            n = 1
            while "{}{}".format ( element, n )[:4] in used or "{}{}".format ( element, n ) in counts: n += 1
            new = "{}{}".format ( element, n )[:4]
        used.add ( new ); names.append ( new )
    connectivity = Connectivity ( )
    for k, ( name, element, xyz, _ ) in enumerate ( atoms ):
        if element not in L.ELEMENT_Z:
            raise ValueError ( "{}: unknown element {}.".format ( path, element ) )
        connectivity.AddNode ( Atom.WithOptions ( atomicNumber = L.ELEMENT_Z[element], label = names[k], formalCharge = charges[k] ) )
    types = { 1: BondType.Single, 2: BondType.Double, 3: BondType.Triple }
    for ( i, j, order ) in kekule:
        connectivity.AddEdge ( Bond.WithNodes ( connectivity.nodes[i], connectivity.nodes[j],
                                                isAromatic = ( min ( i, j ), max ( i, j ) ) in aromatic,
                                                type = types.get ( order, BondType.Single ) ) )
    ConvertInputConnectivity ( connectivity, { } )
    system = System.FromConnectivity ( connectivity = connectivity )
    system.label = label or os.path.splitext ( os.path.basename ( path ) )[0]
    coordinates = Coordinates3.WithExtent ( len ( atoms ) )
    for k, ( _, _, xyz, _ ) in enumerate ( atoms ):
        coordinates[k, 0], coordinates[k, 1], coordinates[k, 2] = float ( xyz[0] ), float ( xyz[1] ), float ( xyz[2] )
    system.coordinates3 = coordinates
    return system


def import_system ( path, log = None ):
    """ ImportSystem(), with system_from_mol2() as the fallback when
        pDynamo's own reader rejects a .mol2 file. Returns ( system,
        used_fallback ). """
    from pBabel import ImportSystem
    try:
        return ImportSystem ( path, log = log ), False
    except Exception as error:
        if not path.lower ( ).endswith ( ".mol2" ) or "connectivity" not in str ( error ).lower ( ):
            raise
        return system_from_mol2 ( path ), True
