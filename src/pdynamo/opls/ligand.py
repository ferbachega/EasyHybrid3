#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: simplified OPLS parametrization of small molecules (phase 4)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Parametrizes a residue that has no PDB component (ligand, cofactor,
#      metal ion) so that prep.py can build and type it:
#
#        1. chemistry: hydrogens at a given pH and bond orders (now in
#           pdynamo/ligand_chemistry.py, force-field-agnostic -- split out
#           2026-11, see that module's header);
#        2. OPLS-AA atom classes by rules (element, hybridization,
#           aromaticity, ring size, neighbours) -- the class names of the
#           OPLS-AA files (CT, CA, C, CM, CZ, CW, CS, CR, NA, NB, NC, N,
#           NT, N2, N3, NO, NZ, O, O2, OH, OS, ON, OY, S, SH, SZ, SY, P,
#           F, Cl, Br, I, HC, HA, H, H3, HO, HS);
#        3. charges: CM5 from GFN1-xTB (the "Mulliken/CM5 charges" block),
#           unscaled by default (CM5_SCALE_NEUTRAL; the Jorgensen group's
#           1.20/1.27 factors were derived for DFT CM5 and over-polarize
#           xTB CM5, see CM5_SCALE_NEUTRAL), averaged over
#           topologically equivalent atoms and corrected to the exact
#           total charge;
#        4. Lennard-Jones per class: the most common (sigma, epsilon) of
#           that class among the OPLS-AA 2008 atom types;
#        5. bonded terms: from the class database (parameters.py) when
#           available, otherwise ESTIMATED -- equilibrium values from the
#           geometry, force constants = median of the database terms with
#           the same elements, missing torsions = 0 -- and listed in the
#           report.
#
#      Output: a LigandParameters object with what prep.py needs: a PDB
#      component file (atoms, bonds, formal charges) for building, and an
#      MM sequence component (atom name -> type, charge) + LJ + type rows
#      for the parameter set. Monoatomic metal ions get their OPLS-AA ion
#      parameters directly (no QC).
#
#      This is a SIMPLIFIED parametrization: types by rules, no torsion
#      fitting, xTB-level CM5 (CM5 was defined on DFT densities). Check the
#      report before production simulations.
#

import os
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict

import numpy as np

from pdynamo.ligand_chemistry import ELEMENT_Z, Z_ELEMENT, OPLS_IONS, LigandParameters , \
                                     find_xtb, find_obabel, _neighbors, _rings, equivalence_classes , \
                                     ligand_chemistry

# . OPLS-AA X-CA-CA-X ring torsion, V2 = 7.25 kcal/mol -> pDynamo k = V2 / 2.
AROMATIC_RING_TORSION_K = 3.625

# . Scale of the GFN1-xTB CM5 charges. 1.20 (the Jorgensen group's universal
#   factor) was derived for CM5 from DFT densities; xTB CM5 is already more
#   polarized: against the OPLS-AA 2008 reference charges of 14 small
#   molecules (tests/benchmarks/opls_ligand_charges.py) the MAE is 0.054 e
#   unscaled, 0.082 e with 1.20, and the least-squares optimum is 0.835.
CM5_SCALE_NEUTRAL = 1.00
CM5_SCALE_CHARGED = 1.00


#=============================================================================
#  OPLS classes
#=============================================================================
def assign_classes ( lig ):
    n = len ( lig.atoms )
    el = [ Z_ELEMENT.get ( z, "X" ) for ( _, z, _ ) in lig.atoms ]
    adj = _neighbors ( n, lig.bonds )
    rings = _rings ( n, adj )
    arom_atoms = set ( )
    for ( i, j ) in lig.aromatic_bonds: arom_atoms.update ( ( i, j ) )
    def ring_size ( i ):
        sizes = [ len ( r ) for r in rings if i in r and r <= arom_atoms ]
        return min ( sizes ) if sizes else None
    def heavy_nbrs ( i ): return [ j for ( j, o ) in adj[i] if el[j] != "H" ]
    def n_h ( i ): return sum ( 1 for ( j, o ) in adj[i] if el[j] == "H" )
    def has_double_to ( i, element ): return any ( o == 2 and el[j] == element for ( j, o ) in adj[i] )
    def max_order ( i ): return max ( [ o for ( _, o ) in adj[i] ] or [ 0 ] )
    classes = [ None ] * n
    q = lig.formal_charges
    # . Heavy atoms.
    for i in range ( n ):
        e = el[i]
        if e == "C":
            if max_order ( i ) == 3 or sum ( 1 for ( _, o ) in adj[i] if o == 2 ) == 2: classes[i] = "CZ"
            elif i in arom_atoms:
                size = ring_size ( i )
                hetero = [ j for j in heavy_nbrs ( i ) if el[j] in ( "N", "O", "S" ) and j in arom_atoms ]
                n_aromatic_rings = sum ( 1 for r in rings if i in r and r <= arom_atoms )
                if n_aromatic_rings >= 2 and any ( len ( r ) == 5 for r in rings if i in r and r <= arom_atoms ):
                    classes[i] = "CB"                       # 5/6 fusion carbon (purine C4/C5, indole C3a/C7a)
                elif size == 5:
                    ring5 = [ r for r in rings if i in r and len ( r ) == 5 and r <= arom_atoms ]
                    fused = any ( len ( r & other ) >= 2 for r in ring5 for other in rings
                                  if other is not r and other <= arom_atoms and other != r )
                    if len ( hetero ) >= 2: classes[i] = "CK" if fused else "CR"   # purine C8 / imidazole C2
                    elif hetero:
                        # . next to a pyridine-like N (2 connections, no H): CV (imidazole/oxazole/
                        #   thiazole C4); next to N-H, O or S: CW (imidazole C5, pyrrole/furan C2)
                        j = hetero[0]
                        classes[i] = "CV" if el[j] == "N" and len ( adj[j] ) == 2 else "CW"
                    else:                   classes[i] = "CS"
                elif len ( [ j for j in hetero if el[j] == "N" ] ) >= 2:
                    classes[i] = "CQ"                       # C between two ring N (pyrimidine/purine C2)
                else:
                    classes[i] = "CA"
            elif has_double_to ( i, "O" ) or has_double_to ( i, "S" ): classes[i] = "C"
            elif any ( o == 2 for ( _, o ) in adj[i] ):
                # . C=N of guanidinium/amidinium: OPLS uses CA
                if has_double_to ( i, "N" ) and sum ( 1 for j in heavy_nbrs ( i ) if el[j] == "N" ) >= 2: classes[i] = "CA"
                else: classes[i] = "CM"
            else: classes[i] = "CT"
        elif e == "N":
            hn = heavy_nbrs ( i )
            if max_order ( i ) == 3: classes[i] = "NZ"
            elif sum ( 1 for j in hn if el[j] == "O" ) >= 2: classes[i] = "NO"
            elif i in arom_atoms:
                ring_partners = { j for r in rings if i in r for j in r }
                if len ( adj[i] ) == 3 and any ( el[j] == "C" and j not in ring_partners and j not in arom_atoms for j in hn ):
                    classes[i] = "N*"                       # aromatic N bonded to an sp3 C (nucleoside N9/N1)
                elif len ( adj[i] ) == 3: classes[i] = "NA"
                else: classes[i] = "NB" if ring_size ( i ) == 5 else "NC"
            elif q[i] == 1 and len ( adj[i] ) == 4: classes[i] = "N3"
            elif any ( el[j] == "C" and any ( o == 2 and el[k] in ( "O", "S" ) for ( k, o ) in adj[j] ) for j in hn ):
                classes[i] = "N"
            elif any ( el[j] == "C" and sum ( 1 for ( k, o ) in adj[j] if el[k] == "N" ) >= 2 and any ( o2 == 2 for ( _, o2 ) in adj[j] ) for j in hn ):
                classes[i] = "N2"
            elif max_order ( i ) == 2: classes[i] = "NB"
            else: classes[i] = "NT"
        elif e == "O":
            hn = heavy_nbrs ( i )
            partner = hn[0] if hn else None
            if n_h ( i ) == 2: classes[i] = "OW"
            elif n_h ( i ) == 1: classes[i] = "OH"
            elif partner is not None and el[partner] == "N" and classes[partner] in ( None, "NO" ) and \
                 sum ( 1 for j in heavy_nbrs ( partner ) if el[j] == "O" ) >= 2: classes[i] = "ON"
            elif partner is not None and el[partner] == "S" and len ( hn ) == 1: classes[i] = "OY"
            elif partner is not None and el[partner] == "P" and len ( hn ) == 1: classes[i] = "O2"
            elif len ( hn ) == 1 and partner is not None and el[partner] == "C":
                # . carboxylate: C bonded to two terminal O, one of them charged
                terminal_o = [ j for ( j, o ) in adj[partner] if el[j] == "O" and len ( adj[j] ) == 1 ]
                if len ( terminal_o ) == 2 and any ( q[j] == -1 for j in terminal_o ): classes[i] = "O2"
                else: classes[i] = "O"
            else: classes[i] = "OS"
        elif e == "S":
            n_o = sum ( 1 for j in heavy_nbrs ( i ) if el[j] == "O" and len ( adj[j] ) == 1 )
            if   n_h ( i ) >= 1: classes[i] = "SH"
            elif n_o >= 2:       classes[i] = "SY"
            elif n_o == 1:       classes[i] = "SZ"
            else:                classes[i] = "S"
        elif e == "P":  classes[i] = "P"
        elif e in ( "F", "Cl", "Br", "I", "Si", "B" ): classes[i] = e
        elif e in OPLS_IONS: classes[i] = e
    # . Hydrogens.
    for i in range ( n ):
        if el[i] != "H": continue
        parent = adj[i][0][0] if adj[i] else None
        pe, pc = ( el[parent], classes[parent] ) if parent is not None else ( None, None )
        if   pe == "C": classes[i] = "HA" if pc in ( "CA", "CW", "CS", "CR", "CV", "CQ", "CB", "CK" ) else "HC"
        elif pe == "N": classes[i] = "H3" if pc == "N3" else "H"
        elif pe == "O": classes[i] = "HW" if pc == "OW" else "HO"
        elif pe == "S": classes[i] = "HS"
        else:           classes[i] = "HC"
    if any ( c is None for c in classes ):
        raise ValueError ( "No OPLS class for atoms {} of {}.".format (
                [ lig.atoms[i][0] for i, c in enumerate ( classes ) if c is None ], lig.residue_name ) )
    lig.classes = classes
    lig.equivalence = equivalence_classes ( el, adj )
    return classes


#=============================================================================
#  Lennard-Jones per class
#=============================================================================
_CLASS_LJ = None

def class_lennard_jones ( ):
    """ { class: ( epsilon, sigma ) } -- most common pair among the
        OPLS-AA 2008 types of the class (united-atom and GBSA types out). """
    global _CLASS_LJ
    if _CLASS_LJ is None:
        from pdynamo.opls.parameters import TINKER_FILES
        from pdynamo.opls.tinker_prm import TinkerPRM
        path = dict ( TINKER_FILES )["oplsaa08"]
        prm = TinkerPRM ( path )
        counts = defaultdict ( Counter )
        for t, data in prm.atoms.items ( ):
            if "(UA)" in data["description"] or "GBSA" in data["description"]: continue
            lj = prm.type_lj ( t )
            if lj is None: continue
            counts[data["symbol"]][( round ( lj[1], 5 ), round ( lj[0], 5 ) )] += 1
        _CLASS_LJ = { s: c.most_common ( 1 )[0][0] for s, c in counts.items ( ) }
        for element, ( _, atom_type ) in OPLS_IONS.items ( ):
            lj = prm.type_lj ( atom_type )
            if lj is not None: _CLASS_LJ["ION_" + element] = ( lj[1], lj[0] )
    return _CLASS_LJ


#=============================================================================
#  Charges
#=============================================================================
def cm5_charges ( lig, xtb = None, scale = None, keep_scratch = None ):
    """ GFN1-xTB CM5 charges (scaled, symmetrized, exact total). """
    xtb = xtb or find_xtb ( )
    if xtb is None:
        raise RuntimeError ( "xtb was not found (set XTBHOME or put xtb in the PATH); it is needed for CM5 charges." )
    total = lig.total_charge
    check_multiplicity ( lig )
    if scale is None: scale = CM5_SCALE_NEUTRAL if total == 0 else CM5_SCALE_CHARGED
    scratch = tempfile.mkdtemp ( prefix = "opls_cm5_" )
    try:
        xyz = os.path.join ( scratch, "lig.xyz" )
        with open ( xyz, "w" ) as handle:
            handle.write ( "{}\n{}\n".format ( len ( lig.atoms ), lig.residue_name ) )
            for ( name, z, r ) in lig.atoms:
                handle.write ( "{} {:.6f} {:.6f} {:.6f}\n".format ( Z_ELEMENT[z], r[0], r[1], r[2] ) )
        command = [ xtb, xyz, "--gfn", "1", "--chrg", str ( total ) ]
        if lig.multiplicity > 1: command += [ "--uhf", str ( lig.multiplicity - 1 ) ]
        result = subprocess.run ( command, cwd = scratch,
                                  capture_output = True, text = True, timeout = 1800,
                                  env = dict ( os.environ, OMP_NUM_THREADS = os.environ.get ( "OMP_NUM_THREADS", "2" ) ) )
        text = result.stdout
        start = text.find ( "Mulliken/CM5 charges" )
        if start < 0:
            raise RuntimeError ( "xtb gave no CM5 charges for {}:\n{}".format ( lig.residue_name, ( result.stdout + result.stderr )[-1500:] ) )
        charges = [ ]
        for line in text[start:].splitlines ( )[1:len ( lig.atoms ) + 1]:
            charges.append ( float ( line.split ( )[2] ) )
        if len ( charges ) != len ( lig.atoms ):
            raise RuntimeError ( "Could not read the CM5 charges of {}.".format ( lig.residue_name ) )
        if keep_scratch:
            shutil.copy ( os.path.join ( scratch, "lig.xyz" ), keep_scratch + ".xyz" )
            with open ( keep_scratch + "_xtb.log", "w" ) as handle: handle.write ( text )
    finally:
        shutil.rmtree ( scratch, ignore_errors = True )
    charges = np.array ( charges ) * scale
    # . Symmetrize over topologically equivalent atoms.
    groups = defaultdict ( list )
    for i, g in enumerate ( lig.equivalence ): groups[g].append ( i )
    for members in groups.values ( ):
        charges[members] = charges[members].mean ( )
    # . Exact total: spread the residual over all atoms.
    charges += ( total - charges.sum ( ) ) / len ( charges )
    # . 6 decimals (as written to the parameter files) with an exact total: the
    #   rounding residual (< 1e-5) goes to the atom with the largest |q|
    charges = np.round ( charges, 6 )
    charges[int ( np.argmax ( np.abs ( charges ) ) )] += round ( total - float ( charges.sum ( ) ), 6 )
    lig.charges = [ float ( c ) for c in charges ]
    lig.charge_method = "CM5 (GFN1-xTB) x {:.2f}".format ( scale )
    if lig.multiplicity > 1:
        lig.charge_method += ", multiplicity {}".format ( lig.multiplicity )
    return lig.charges


#=============================================================================
#  Labels, LJ
#=============================================================================
def finalize_types ( lig ):
    """ Type labels CLASS.RESn (one per equivalence class) and LJ. """
    lj = class_lennard_jones ( )
    if len ( lig.atoms ) == 1 and lig.notes and lig.notes[0].startswith ( "monoatomic" ):
        element = Z_ELEMENT[lig.atoms[0][1]]
        label = "{}.{}".format ( element, lig.residue_name )
        lig.classes, lig.labels = [ element ], [ label ]
        lig.charges = [ float ( lig.formal_charges[0] ) ]
        lig.lj = { label: lj["ION_" + element] }
        lig.charge_method = "formal charge (OPLS-AA ion)"
        lig.equivalence = [ 0 ]
        return
    labels = [ ]
    for i, ( c, g ) in enumerate ( zip ( lig.classes, lig.equivalence ) ):
        labels.append ( "{}.{}{}".format ( c, lig.residue_name, g + 1 ) )
    lig.labels = labels
    for label, c in zip ( labels, lig.classes ):
        if c not in lj:
            lig.warnings.append ( "no LJ for class {} -- using CT values".format ( c ) )
        lig.lj[label] = lj.get ( c, lj["CT"] )

def check_multiplicity ( lig ):
    """ Raises ValueError when the number of electrons (sum of Z minus the
        total charge) and the multiplicity are incompatible (an even number
        of electrons needs an odd multiplicity and vice versa). """
    m = int ( lig.multiplicity )
    if m < 1:
        raise ValueError ( "{}: the multiplicity must be 1 or larger.".format ( lig.residue_name ) )
    if ( lig.n_electrons + m - 1 ) % 2 != 0:
        raise ValueError ( "{}: charge {:+d} gives {} electrons, incompatible with multiplicity {} "
                           "({} electrons need an {} multiplicity).".format (
                           lig.residue_name, lig.total_charge, lig.n_electrons, m, "even" if lig.n_electrons % 2 == 0 else "odd",
                           "odd" if lig.n_electrons % 2 == 0 else "even" ) )

def apply_user_charge ( lig, total_charge = None, multiplicity = None ):
    """ Total charge / multiplicity chosen by the user. The formal charges
        of the atoms (used for building and for the atom classes) are not
        changed: a different total only changes the number of electrons of
        the QC charge calculation and the normalization of the charges. To
        change protons, change the pH or the hydrogens instead. """
    if total_charge is not None and int ( total_charge ) != lig.formal_charge_sum:
        lig.charge_override = int ( total_charge )
        lig.warnings.append ( "total charge {:+d} set by the user (perceived formal charges sum to {:+d})".format (
                lig.charge_override, lig.formal_charge_sum ) )
    if multiplicity is not None and int ( multiplicity ) != 1:
        lig.multiplicity = int ( multiplicity )
        lig.notes.append ( "spin multiplicity {} set by the user".format ( lig.multiplicity ) )
    check_multiplicity ( lig )
    return lig

def parametrize_residue ( residue_name, atoms, ph = 7.4, xtb = None, cm5_scale = None, add_hydrogens = None,
                          keep_scratch = None, bonds = None, formal_charges = None, total_charge = None,
                          multiplicity = None, h_counts = None ):
    """ Full simplified parametrization of one residue. total_charge /
        multiplicity: the user's values (see apply_user_charge()). """
    lig = ligand_chemistry ( residue_name, atoms, ph = ph, add_hydrogens = add_hydrogens, bonds = bonds,
                             formal_charges = formal_charges, h_counts = h_counts )
    if len ( lig.atoms ) == 1 and ( total_charge is not None or multiplicity not in ( None, 1 ) ):
        lig.warnings.append ( "charge/multiplicity of a monoatomic ion come from the OPLS-AA ion parameters (user values ignored)" )
    elif len ( lig.atoms ) > 1:
        apply_user_charge ( lig, total_charge, multiplicity )
    if len ( lig.atoms ) > 1:
        assign_classes ( lig )
        cm5_charges ( lig, xtb = xtb, scale = cm5_scale, keep_scratch = keep_scratch )
    finalize_types ( lig )
    return lig


#=============================================================================
#  Missing bonded terms
#=============================================================================
def estimate_missing_terms ( system, types, database, base_set = None ):
    """ { term: { key: value } } for the bonds/angles/dihedrals of `system`
        that neither the base set nor the class database has: equilibrium
        values from the current geometry, force constants = median over the
        database terms with the same element tuple (fallback 300 / 50),
        torsions = no term (empty list). Also returns notes. """
    from pdynamo.opls import parameters as P
    used = P.system_terms ( system, types )
    crd  = system.coordinates3
    z    = [ a.atomicNumber for a in system.atoms ]
    element_of_class = { }
    for label, z_ in zip ( types, z ): element_of_class[P.atom_class ( label )] = z_
    def elements ( key ): return tuple ( sorted ( element_of_class.get ( P.atom_class ( l ), 0 ) for l in key ) )
    medians = { "bond": defaultdict ( list ), "angle": defaultdict ( list ) }
    for term in ( "bond", "angle" ):
        for key, entry in database.entries[term].items ( ):
            medians[term][elements ( key )].append ( entry["value"][1] )
    # . geometry per key
    geometry = { "bond": defaultdict ( list ), "angle": defaultdict ( list ) }
    c = system.connectivity
    def xyz ( i ): return np.array ( [ crd[i, 0], crd[i, 1], crd[i, 2] ] )
    for ( i, j ) in c.bondIndices:
        geometry["bond"][P.key_bond ( types[i], types[j] )].append ( float ( np.linalg.norm ( xyz ( i ) - xyz ( j ) ) ) )
    for ( i, j, k ) in c.angleIndices:
        u, v = xyz ( i ) - xyz ( j ), xyz ( k ) - xyz ( j )
        cosine = np.dot ( u, v ) / ( np.linalg.norm ( u ) * np.linalg.norm ( v ) )
        geometry["angle"][P.key_angle ( types[i], types[j], types[k] )].append ( float ( np.degrees ( np.arccos ( np.clip ( cosine, -1, 1 ) ) ) ) )
    # . dihedrals whose central bond is aromatic (both central atoms aromatic, in one ring)
    aromatic_central = set ( )
    for ( i, j, k, l ) in c.dihedralIndices:
        if getattr ( system.atoms[j], "isAromatic", False ) and getattr ( system.atoms[k], "isAromatic", False ):
            aromatic_central.add ( P.key_dihedral ( types[i], types[j], types[k], types[l] ) )
    out = { "bond": { }, "angle": { }, "dihedral": { }, "outofplane": { } }
    notes = [ ]
    # . out-of-plane terms: OPLS-AA impropers for planar centres (Tinker V / 2):
    #   carbonyl C 10.5, aromatic C/N 1.1, amide-like N 1.0; others stay 0
    central_aromatic = { }
    central_carbonyl = { }
    for ( i, j, k, l ) in c._IndicesNeighbor3 ( ):
        key = P.key_outofplane ( types[i], types[j], types[k], types[l] )
        central_aromatic[key] = bool ( getattr ( system.atoms[i], "isAromatic", False ) )
        central_carbonyl[key] = system.atoms[i].atomicNumber == 6 and any (
                system.atoms[n].atomicNumber == 8 and len ( c.adjacentNodes[system.atoms[n]] ) == 1 for n in ( j, k, l ) )
    for key in sorted ( used["outofplane"] ):
        if base_set is not None and key in base_set.terms["outofplane"]: continue
        value, _, _ = database.lookup ( "outofplane", key )
        if value is not None: continue
        central = P.atom_class ( key[0] )
        if central_carbonyl.get ( key ):      K = 10.5
        elif central_aromatic.get ( key ):    K = 1.1
        elif central in ( "N", "N2", "NA", "N*" ): K = 1.0
        else: continue
        out["outofplane"][key] = K
        notes.append ( "outofplane {} : generic OPLS improper (K {:.1f})".format ( "-".join ( key ), K ) )
    for term in ( "bond", "angle", "dihedral" ):
        for key in sorted ( used[term] ):
            if base_set is not None and key in base_set.terms[term]: continue
            value, _, _ = database.lookup ( term, key )
            if value is not None: continue
            if term == "dihedral":
                if key in aromatic_central:
                    out[term][key] = [ ( AROMATIC_RING_TORSION_K, 2, 180.0 ) ]
                    notes.append ( "dihedral {} : generic aromatic ring torsion (k {:.3f}, n 2)".format ( "-".join ( key ), AROMATIC_RING_TORSION_K ) )
                else:
                    out[term][key] = [ ]
                    notes.append ( "dihedral {} : no parameters, set to zero".format ( "-".join ( key ) ) )
                continue
            ks = medians[term].get ( elements ( key ) )
            k_value = float ( np.median ( ks ) ) if ks else ( 300.0 if term == "bond" else 50.0 )
            eq = float ( np.mean ( geometry[term][key] ) )
            out[term][key] = ( round ( eq, 4 ), round ( k_value, 2 ) )
            notes.append ( "{} {} : estimated ({:.3f}, k {:.1f})".format ( term, "-".join ( key ), eq, k_value ) )
    return out, notes
