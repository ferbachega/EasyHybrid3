#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: solvation of a READY system (one that already has its MM model)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description (2026-10-21, user's request: "uma ferramenta para solvatar
#  sistemas depois de prontos", replacing Edit > Solvate System, which needed a
#  pre-built solvent-box system and ion .pkl files):
#
#    1. Placement: the Builder's solvation core (gui/windows/builder/
#       solvation.py: replicated pre-equilibrated TIP3P box, seam clashes
#       removed, waters too close to the solute removed, counter-ions and salt
#       placed by replacing waters).
#    2. The waters and ions become a pDynamo system of their own (PDB route
#       with the component library: HOH, NA, K, CL... with bonds and formal
#       charges) and get the SAME force field as the solute:
#         OPLS  the solute's own parameter set (TIP3P water, OPLS ions);
#         DYFF  the solute's DYFF set, TIP3P charges, the TIP3P oxygen LJ and
#               the polar-hydrogen LJ type, ions at their formal charge --
#               what the residue library gives water and ions
#               (residue_library.solvent_and_ions, dyff_edits.water_lj).
#    3. pDynamo's MergeByAtom joins solute (first: its atom indices do not
#       change) and solvent; their MM states (terms, charges, LJ types) are
#       merged as they are, so the solute keeps exactly its parameters.
#    4. A periodic box gets cubic/orthorhombic symmetry and a cut-off NB
#       model; a droplet stays non-periodic.
#
#  AMBER and CHARMM systems are not handled: their water and ions must come
#  from the same topology program (tleap, CHARMM-GUI/psfgen).
#

import os
import tempfile

import numpy as np

from pCore import Clone

SUPPORTED = { "OPLS": "OPLS", "DYFF": "DYFF" }
ION_RESIDUES = { "Na+": ( "NA", "NA" ), "K+": ( "K", "K" ), "Li+": ( "LI", "LI" ), "Cl-": ( "CL", "CL" ), "Br-": ( "BR", "BR" ) }
IONS_BY_FORCE_FIELD = { "OPLS": ( ( "Na+", "K+" ), ( "Cl-", ) ), "DYFF": ( ( "Na+", "K+", "Li+" ), ( "Cl-", "Br-" ) ) }
WATERS_PER_CHAIN = 9999
_WATER_CHAINS = "WXYZUVTSRQPONMLKJHGFEDCBA"


class SolvateSystemError ( Exception ):
    pass


def force_field_kind ( system ):
    """ "OPLS", "DYFF" or None. """
    model = getattr ( system, "mmModel", None )
    if model is None: return None
    return SUPPORTED.get ( str ( getattr ( model, "forceField", "" ) ).upper ( ) )


def check ( system ):
    """ Raises SolvateSystemError with a message for the user when the system cannot be solvated here. """
    if system is None:
        raise SolvateSystemError ( "No system selected." )
    model = getattr ( system, "mmModel", None )
    if model is None:
        raise SolvateSystemError ( "'{}' has no force field. Define one first (System > Force Field > Define Force Field): "
                                   "the water and the ions are given the solute's force field.".format ( system.label ) )
    if force_field_kind ( system ) is None:
        raise SolvateSystemError ( "'{}' uses {}: solvate it with the program that wrote its topology (tleap for AMBER, "
                                   "CHARMM-GUI or psfgen for CHARMM). This tool solvates OPLS and DYFF systems.".format (
                                   system.label, getattr ( model, "forceField", type ( model ).__name__ ) ) )
    if getattr ( system, "qcModel", None ) is not None:
        raise SolvateSystemError ( "'{}' has a QC region. Solvate first, then define the QC region on the solvated "
                                   "system (System > Remove > QC Atoms to remove it).".format ( system.label ) )
    if getattr ( system, "symmetry", None ) is not None:
        raise SolvateSystemError ( "'{}' already has a periodic box (is it already solvated?).".format ( system.label ) )


def solute_charge ( system ):
    """ ( rounded total charge, exact total ) from the MM charges. """
    total = float ( sum ( system.mmState.charges ) )
    return int ( round ( total ) ), total


def system_xyz ( system, coordinates3 = None ):
    crd = coordinates3 if coordinates3 is not None else system.coordinates3
    return np.array ( [ [ crd[i, 0], crd[i, 1], crd[i, 2] ] for i in range ( len ( system.atoms ) ) ], dtype = float )


def water_template ( ):
    from gui.windows.builder import solvation
    return solvation.default_solvent ( )


def place ( xyz, charge, shape = "cube", padding = 10.0, box_size = None, neutralize = True, concentration = 0.15,
            cation = "Na+", anion = "Cl-", closeness = 2.4, align = False, seed = 1 ):
    """ The Builder's solvation.solvate() on the solute coordinates (fast, numpy only). """
    from gui.windows.builder import solvation
    try:
        return solvation.solvate ( np.asarray ( xyz, dtype = float ), solute_charge = charge, shape = shape, padding = padding,
                                   box_size = box_size, closeness = closeness, neutralize = neutralize,
                                   concentration = concentration, cation = cation, anion = anion, align = align,
                                   seed = seed, solvent = water_template ( ) )
    except solvation.SolvationError as error:
        raise SolvateSystemError ( str ( error ) )


# ----------------------------------------------------------------------------
#  solvent system
# ----------------------------------------------------------------------------
def _pdb_atom ( serial, name, residue, chain, number, xyz, element ):
    name = name if len ( name ) >= 4 else " " + name.ljust ( 3 )
    return "HETATM{:5d} {:4s} {:>3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00          {:>2s}".format (
           serial % 100000, name, residue, chain, number % 10000, xyz[0], xyz[1], xyz[2], element.upper ( ) )


def write_solvent_pdb ( result, path ):
    template = water_template ( )
    lines, serial = [ ], 0
    for k, molecule in enumerate ( result.waters ):
        chain = _WATER_CHAINS[k // WATERS_PER_CHAIN]
        for name, element, xyz in zip ( template.names, template.elements, molecule ):
            serial += 1
            lines.append ( _pdb_atom ( serial, name, "HOH", chain, k % WATERS_PER_CHAIN + 1, xyz, element ) )
    for k, ( label, element, name, charge, xyz ) in enumerate ( result.ions ):
        serial += 1
        residue, atom_name = ION_RESIDUES[label]
        lines.append ( _pdb_atom ( serial, atom_name, residue, "I", k + 1, xyz, element ) )
    with open ( path, "w" ) as handle:
        handle.write ( "\n".join ( lines ) + "\nEND\n" )


def build_solvent_system ( result, solute ):
    """ The waters and ions of `result` as a pDynamo system with the solute's force field. """
    import contextlib
    import io
    from pBabel import ImportSystem
    if result.n_waters == 0 and not result.ions:
        return None
    fd, path = tempfile.mkstemp ( suffix = ".pdb", prefix = "easyhybrid_solvent_" ); os.close ( fd )
    try:
        write_solvent_pdb ( result, path )
        with contextlib.redirect_stdout ( io.StringIO ( ) ):
            solvent = ImportSystem ( path, useComponentLibrary = True, log = None )
    finally:
        os.remove ( path )
    expected = 3 * result.n_waters + len ( result.ions )
    if len ( solvent.atoms ) != expected:
        raise SolvateSystemError ( "The solvent was read with {} atoms instead of {}.".format ( len ( solvent.atoms ), expected ) )
    kind = force_field_kind ( solute )
    if kind == "OPLS":
        solvent.DefineMMModel ( Clone ( solute.mmModel ), log = None )
    else:
        _dyff_solvent ( solvent, solute )
    # . LJ containers merge only when their labels match
    for name in ( "ljParameters", "ljParameters14" ):
        mine, theirs = getattr ( solvent.mmState, name, None ), getattr ( solute.mmState, name, None )
        if mine is not None and theirs is not None:
            mine.label = theirs.label
    return solvent


def _dyff_solvent ( solvent, solute ):
    from pdynamo.dyff_edits import dyff_model, water_lj, soften_polar_hydrogens
    from pdynamo.residue_library import WATER_CHARGES
    model = solute.mmModel
    name = getattr ( model, "parameterSet", None ) or "dyff-1.0"
    solvent.DefineMMModel ( dyff_model ( name ), log = None )
    charges = solvent.mmState.charges
    for i, atom in enumerate ( solvent.atoms ):
        if atom.atomicNumber == 8:   charges[i] = WATER_CHARGES["O"]
        elif atom.atomicNumber == 1: charges[i] = WATER_CHARGES["H"]
        else:                        charges[i] = float ( atom.formalCharge or 0 )
    water_lj ( solvent )
    soften_polar_hydrogens ( solvent )
    if solvent.mmModel.forceField != model.forceField:
        raise SolvateSystemError ( "The DYFF model of the solvent ({}) does not match the solute's ({}).".format (
                                   solvent.mmModel.forceField, model.forceField ) )


# ----------------------------------------------------------------------------
#  the whole operation
# ----------------------------------------------------------------------------
def apply_cell ( system, result ):
    if result.cell is None: return
    from pScientific.Symmetry import PeriodicBoundaryConditions, CrystalSystemCubic, CrystalSystemOrthorhombic
    a, b, c = result.cell[:3]
    if abs ( a - b ) < 1e-6 and abs ( a - c ) < 1e-6:
        crystal, parameters = CrystalSystemCubic ( ), { "a": a }
    else:
        crystal, parameters = CrystalSystemOrthorhombic ( ), { "a": a, "b": b, "c": c }
    system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( crystal )
    system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( **parameters )


def solvate ( system, coordinates3 = None, label = None, progress = None, **options ):
    """ ( solvated system, report ). options: see place(). The solute is the
        first block of atoms of the new system, with its MM parameters unchanged. """
    import contextlib
    import io
    from pScientific.Geometry3 import Coordinates3
    from pSimulation import MergeByAtom
    from pMolecule.NBModel import NBModelCutOff
    say = progress or ( lambda text: None )
    check ( system )
    kind = force_field_kind ( system )
    cations, anions = IONS_BY_FORCE_FIELD[kind]
    if options.get ( "cation", "Na+" ) not in cations or options.get ( "anion", "Cl-" ) not in anions:
        raise SolvateSystemError ( "With {} the ions available are {} and {}.".format ( kind, ", ".join ( cations ), ", ".join ( anions ) ) )
    charge, exact = solute_charge ( system )
    say ( "placing water and ions" )
    result = place ( system_xyz ( system, coordinates3 ), charge, **options )
    say ( "{} waters, {} ions: assigning {} parameters".format ( result.n_waters, len ( result.ions ), kind ) )
    solvent = build_solvent_system ( result, system )
    # . the solute, moved into the box, without the parts that cannot be merged
    solute = Clone ( system )
    crd = Coordinates3.WithExtent ( len ( solute.atoms ) )
    for i, xyz in enumerate ( result.solute_xyz ): crd[i, 0], crd[i, 1], crd[i, 2] = xyz
    solute.coordinates3 = crd
    nb_model = getattr ( system, "nbModel", None )
    solute.DefineRestraintModel ( None )
    solute.freeAtoms = None
    solute.nbModel = None
    say ( "merging" )
    with contextlib.redirect_stdout ( io.StringIO ( ) ):
        new = MergeByAtom ( solute, solvent ) if solvent is not None else solute
    if new.mmState is None or new.mmState.ljParameters is None:
        raise SolvateSystemError ( "pDynamo could not merge the solute and solvent MM models." )
    new.label = label or "{} (solvated)".format ( system.label )
    apply_cell ( new, result )
    if result.cell is not None:
        keep = nb_model is not None and type ( nb_model ).__name__ in ( "NBModelCutOff", "NBModelABFS" )
        new.DefineNBModel ( Clone ( nb_model ) if keep else NBModelCutOff.WithDefaults ( ) )
    else:
        new.DefineNBModel ( Clone ( nb_model ) if nb_model is not None else NBModelCutOff.WithDefaults ( ) )
    ions = { }
    for ( ion, *_ ) in result.ions: ions[ion] = ions.get ( ion, 0 ) + 1
    final_charge = float ( sum ( new.mmState.charges ) )
    report = { "force_field": kind, "shape": result.region.shape, "cell": result.cell,
               "sphere": ( [ float ( v ) for v in result.sphere[0] ], float ( result.sphere[1] ) ) if result.sphere else None,
               "waters": result.n_waters, "ions": ions, "solute_charge": charge, "solute_charge_exact": exact,
               "final_charge": final_charge, "removed_by_solute": result.removed_by_solute, "atoms": len ( new.atoms ),
               "solute_atoms": len ( system.atoms ), "notes": list ( result.notes ) }
    return new, report


def report_text ( report ):
    lines = [ "Force field: {} (TIP3P water)".format ( report["force_field"] ) ]
    if report["cell"]:
        lines.append ( "Box: {:.2f} x {:.2f} x {:.2f} Å ({})".format ( *report["cell"][:3], report["shape"] ) )
    else:
        lines.append ( "Droplet: radius {:.2f} Å (no periodic box)".format ( report["sphere"][1] ) )
    lines.append ( "Waters: {} ({} removed for overlapping the solute)".format ( report["waters"], report["removed_by_solute"] ) )
    lines.append ( "Ions: {}".format ( ", ".join ( "{} {}".format ( n, k ) for k, n in sorted ( report["ions"].items ( ) ) ) or "none" ) )
    lines.append ( "Charge: solute {:+.3f}, solvated system {:+.3f}".format ( report["solute_charge_exact"], report["final_charge"] ) )
    lines.append ( "Atoms: {} ({} solute + {} solvent)".format ( report["atoms"], report["solute_atoms"], report["atoms"] - report["solute_atoms"] ) )
    lines += report.get ( "notes", [ ] )
    return "\n".join ( lines )
