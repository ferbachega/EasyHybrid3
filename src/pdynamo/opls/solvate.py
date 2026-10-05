#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: OPLS system preparation -- solvation (phase 2, no GUI)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Solvates a PREPARED system (prep.build_prepared_system(), hydrogens
#      included) with the Builder's solvation core (gui/windows/builder/
#      solvation.py: replicated pre-equilibrated TIP3P box, seam clash
#      removal, ion replacement) and turns waters/ions into residues of
#      a new PrepAnalysis, so that the final system is built by the SAME
#      pDynamo PDB-model route as the solute: HOH/NA/K/CL components with
#      their library names, bonds and formal charges, typed by the OPLS
#      "Water"/ion patterns (TIP3P charges, OW -0.834 / HW 0.417).
#
#      Periodic shapes (cube, box) give the system a cubic/orthorhombic
#      symmetry; the sphere is a non-periodic droplet.
#

import numpy as np

from pdynamo.opls import prep

# . Ions the OPLS "protein" parameter set types (patterns + LJ).
OPLS_CATIONS = ( "Na+", "K+" )
OPLS_ANIONS  = ( "Cl-", )
_ION_RESIDUE = { "Na+": ( "NA", "NA", 11 ), "K+": ( "K", "K", 19 ), "Cl-": ( "CL", "CL", 17 ) }

WATER_SEGMENT_SIZE = 9999        # residues per water segment (PDB resSeq has 4 digits)


def _water_template ( ):
    from gui.windows.builder import solvation
    return solvation.default_solvent ( )

def solvate_coordinates ( solute_xyz, solute_charge, shape = "cube", padding = 10.0, box_size = None,
                          neutralize = True, concentration = 0.15, cation = "Na+", anion = "Cl-",
                          closeness = 2.4, align = False, seed = 1 ):
    """ Thin wrapper over the Builder's solvation.solvate() restricted to
        the OPLS-supported ions. Returns its SolvationResult. """
    from gui.windows.builder import solvation
    if cation not in OPLS_CATIONS:
        raise ValueError ( "Cation {} has no OPLS parameters here (use {}).".format ( cation, ", ".join ( OPLS_CATIONS ) ) )
    if anion not in OPLS_ANIONS:
        raise ValueError ( "Anion {} has no OPLS parameters here (use {}).".format ( anion, ", ".join ( OPLS_ANIONS ) ) )
    return solvation.solvate ( np.asarray ( solute_xyz, dtype = float ), solute_charge = solute_charge,
                               shape = shape, padding = padding, box_size = box_size, closeness = closeness,
                               neutralize = neutralize, concentration = concentration, cation = cation,
                               anion = anion, align = align, seed = seed, solvent = _water_template ( ) )

def _system_xyz ( system ):
    crd = system.coordinates3
    return np.array ( [ [ crd[i, 0], crd[i, 1], crd[i, 2] ] for i in range ( len ( system.atoms ) ) ], dtype = float )

def _coordinates3 ( xyz ):
    from pScientific.Geometry3 import Coordinates3
    crd = Coordinates3.WithExtent ( len ( xyz ) )
    for i, ( x, y, z ) in enumerate ( xyz ):
        crd[i, 0], crd[i, 1], crd[i, 2] = x, y, z
    return crd

def _residue ( name, number, segment, atoms ):
    r = prep.PrepResidue ( )
    r.entity, r.name, r.number, r.library_name = segment, name, number, name
    r.kind = "water" if name == "HOH" else "ion"
    r.atoms = [ ( atom_name, z, np.asarray ( xyz, dtype = float ) ) for ( atom_name, z, xyz ) in atoms ]
    r.n_hydrogens = sum ( 1 for a in r.atoms if a[1] == 1 )
    r.segment = segment
    return r

def _free_segment ( analysis, base ):
    used = set ( analysis.segments.keys ( ) )
    for k in range ( 1, 10000 ):
        label = "{}{}".format ( base, k )[:4]
        if label not in used: return label
    raise ValueError ( "No free segment label." )

def solvated_analysis ( prepared_system, result, ligands = None, extra_library_paths = None ):
    """ A PrepAnalysis of the prepared solute (at its solvated position)
        plus one residue per water/ion of `result`. """
    analysis = prep.analyze_system ( prepared_system, coordinates3 = _coordinates3 ( result.solute_xyz ),
                                     known_ligands = ligands, extra_library_paths = extra_library_paths )
    template = _water_template ( )
    elements = { "O": 8, "H": 1 }
    # . Waters, in segments of at most WATER_SEGMENT_SIZE residues.
    segment, count = None, WATER_SEGMENT_SIZE
    for molecule in result.molecules:
        if count >= WATER_SEGMENT_SIZE:
            segment, count = _free_segment ( analysis, "W" ), 0
            analysis.segments[segment] = [ ]
        count += 1
        atoms = [ ( name, elements[symbol], xyz ) for name, symbol, xyz in zip ( template.names, template.elements, molecule ) ]
        analysis.residues.append ( _residue ( "HOH", count, segment, atoms ) )
        analysis.segments[segment].append ( len ( analysis.residues ) - 1 )
    # . Ions.
    if result.ions:
        segment = _free_segment ( analysis, "I" )
        analysis.segments[segment] = [ ]
        for k, ( label, element, name, charge, xyz ) in enumerate ( result.ions ):
            residue_name, atom_name, z = _ION_RESIDUE[label]
            analysis.residues.append ( _residue ( residue_name, k + 1, segment, [ ( atom_name, z, xyz ) ] ) )
            analysis.segments[segment].append ( len ( analysis.residues ) - 1 )
    for r in analysis.residues:
        r.missing_heavy = [ ]
    return analysis

def apply_cell ( system, result ):
    """ Cubic/orthorhombic symmetry from a periodic SolvationResult. """
    if result.cell is None: return
    from pScientific.Symmetry import PeriodicBoundaryConditions, CrystalSystemCubic, CrystalSystemOrthorhombic
    a, b, c = result.cell[:3]
    if abs ( a - b ) < 1e-6 and abs ( a - c ) < 1e-6:
        crystal, parameters = CrystalSystemCubic ( ), { "a": a }
    else:
        crystal, parameters = CrystalSystemOrthorhombic ( ), { "a": a, "b": b, "c": c }
    system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( crystal )
    system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( **parameters )

def solvate_prepared_system ( prepared_system, label = None, ligands = None, extra_library_paths = None,
                              charge_correction = 0, **options ):
    """ ( solvated system without MM model, report ). options: see
        solvate_coordinates(). The solute charge is the sum of the
        prepared system's formal charges plus charge_correction (ligand
        total charges set by the user, see ligand.apply_user_charge()). """
    charge = int ( sum ( getattr ( a, "formalCharge", 0 ) or 0 for a in prepared_system.atoms ) ) + int ( charge_correction )
    result = solvate_coordinates ( _system_xyz ( prepared_system ), charge, **options )
    analysis = solvated_analysis ( prepared_system, result, ligands = ligands, extra_library_paths = extra_library_paths )
    system, report = prep.build_prepared_system ( analysis, label = label or getattr ( prepared_system, "label", None ) )
    apply_cell ( system, result )
    ions = { }
    for ( ion_label, *_ ) in result.ions:
        ions[ion_label] = ions.get ( ion_label, 0 ) + 1
    report.update ( { "shape": result.region.shape, "cell": result.cell,
                      "sphere": ( [ float ( v ) for v in result.sphere[0] ], float ( result.sphere[1] ) ) if result.sphere else None,
                      "waters": result.n_waters, "ions": ions, "solute_charge": charge,
                      "removed_by_solute": result.removed_by_solute } )
    return system, report
