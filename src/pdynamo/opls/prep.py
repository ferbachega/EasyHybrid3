#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: OPLS system preparation -- phase 0 (no GUI)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      OPLS-specific tail of what used to be the whole of this file (split
#      2026-11: the generic residue/ligand/system analysis and building
#      engine -- analyze_system, build_prepared_system, write_clean_pdb,
#      library_paths, ... -- moved to pdynamo/structure_prep.py, since
#      DYFF's own preparation needed it too and none of it types anything
#      with OPLS; see that module's header for the full rationale). What is
#      left here is specific to the OPLS-AA force field itself:
#
#        type_atoms ( system )               -> ( types, charges, untyped )
#        typing_report ( system )             -> dict
#        prepare_opls_system ( system )       -> ( new System, report )
#        parametrize_ligands ( analysis, ... ) -> assigns OPLS-AA atom
#            classes (pdynamo/opls/ligand.py) to every unsupported residue,
#            unlike structure_prep.ligand_preview/ligand_protonation_baseline,
#            which only need the force-field-agnostic chemistry step.
#
#      Residues without a PDB library component (ligands, cofactors) are
#      parametrized separately (parametrize_ligands) before build.
#

import os
import tempfile
from collections import Counter, OrderedDict

import numpy as np

from pdynamo.structure_prep import ( analyze_system, build_prepared_system, ligand_preview , \
                                     ligand_protonation_baseline, ligands_to_parametrize, library_paths , \
                                     PrepAnalysis, PrepResidue, PROTONATION_CHOICES, source_atom_map , \
                                     split_ligand_variants, truncate_residue, write_clean_pdb )

# . Default OPLS parameter set (pDynamo3's own).
DEFAULT_PARAMETER_SET = "protein"


#=============================================================================
#  Ligands (phase 4)
#=============================================================================
COVALENT_CONTACT = 2.0      # A, ligand heavy atom closer than this to another residue = covalent


class CloseContactError ( ValueError ):
    """ [EN] 2026-10-17: a residue to parametrize has a heavy atom closer than
        COVALENT_CONTACT to another residue WITHOUT a bond between them -- maybe
        a covalent link missing from the input. The user can ignore it
        (prepare_opls_system ( trust_connectivity = True )): only the bonds of
        the input are bonds. """
    def __init__ ( self, message, residue = "", atom = "", other = "", distance = 0.0 ):
        super ( ).__init__ ( message )
        self.residue, self.atom, self.other, self.distance = residue, atom, other, distance

def parametrize_ligands ( analysis, work_dir, ph = 7.4, xtb = None, cm5_scale = None, progress = None,
                          trust_connectivity = False ):
    """ Parametrizes every unsupported residue of the analysis
        (pdynamo/opls/ligand.py: Open Babel hydrogens at pH, OPLS classes by
        rules, GFN1-xTB CM5 charges) and makes it buildable: its PDB
        component goes to <work_dir>/pdb_components, the residue becomes a
        "ligand". One parametrization per residue NAME; every copy gets its
        own hydrogens. Raises CloseContactError for a heavy atom closer than
        COVALENT_CONTACT to another residue, unless trust_connectivity (the
        input's bonds are the only bonds: a Builder system). """
    from pdynamo.opls import ligand as L
    say = progress or ( lambda text: None )
    if not analysis.unsupported: return { }
    residues = analysis.residues
    heavy_others = [ ( i, xyz, n ) for i, r in enumerate ( residues ) for ( n, z, xyz ) in r.atoms
                     if z != 1 and xyz is not None and r.kind != "water" ]
    others_xyz = np.array ( [ x for ( _, x, _ ) in heavy_others ] ) if heavy_others else np.zeros ( ( 0, 3 ) )
    others_id  = np.array ( [ i for ( i, _, _ ) in heavy_others ] ) if heavy_others else np.zeros ( 0, int )
    others_nm  = [ n for ( _, _, n ) in heavy_others ]
    library_root = os.path.join ( work_dir, "pdb_components" )
    split_ligand_variants ( analysis )
    by_name = OrderedDict ( )
    for i in analysis.unsupported:
        by_name.setdefault ( residues[i].name, [ ] ).append ( i )
    # . every contact checked before anything is changed: after a CloseContactError
    #   the analysis is intact and the caller can run again with trust_connectivity
    for name, indices in ( [ ] if trust_connectivity else by_name.items ( ) ):
        for i in indices:
            r = residues[i]
            mine_atoms = [ ( n, xyz ) for ( n, z, xyz ) in r.atoms if z != 1 and xyz is not None ]
            mine = np.array ( [ xyz for ( n, xyz ) in mine_atoms ] )
            if len ( mine ) and len ( others_xyz ):
                mask = others_id != i
                d = np.linalg.norm ( others_xyz[mask][:, None, :] - mine[None, :, :], axis = 2 )
                if d.size and d.min ( ) < COVALENT_CONTACT and len ( r.atoms ) > 1:
                    k_other, k_mine = np.unravel_index ( int ( np.argmin ( d ) ), d.shape )
                    other_index = np.nonzero ( mask )[0][k_other]
                    other = "{} {}".format ( residues[others_id[other_index]].label, others_nm[other_index] )
                    atom  = "{} {}".format ( r.label, mine_atoms[k_mine][0] )
                    raise CloseContactError ( "{} and {} are {:.2f} A apart (closer than {:.1f} A: covalent "
                                              "link?); covalent ligands are not supported.".format (
                                              atom, other, float ( d.min ( ) ), COVALENT_CONTACT ),
                                              residue = r.label, atom = atom, other = other, distance = float ( d.min ( ) ) )
    for name, indices in by_name.items ( ):
        first = residues[indices[0]]
        say ( "Parametrizing {} (CM5 charges)...".format ( name ) )
        options = analysis.ligand_options.get ( name, { } )
        lig = L.parametrize_residue ( name, first.atoms, ph = ph, xtb = xtb, cm5_scale = cm5_scale,
                                      keep_scratch = os.path.join ( work_dir, "ligand_" + name ),
                                      bonds = first.bonds, formal_charges = first.formal_charges,
                                      total_charge = options.get ( "charge" ), multiplicity = options.get ( "multiplicity" ),
                                      h_counts = options.get ( "h_counts" ) )
        lig.write_pdb_component ( library_root )
        names = [ a[0] for a in lig.atoms ]
        for i in indices:
            r = residues[i]
            atoms = lig.atoms if i == indices[0] else L.ligand_chemistry ( name, r.atoms, ph = ph, bonds = r.bonds,
                                                                             formal_charges = r.formal_charges,
                                                                             h_counts = options.get ( "h_counts" ) ).atoms
            if [ a[0] for a in atoms ] != names:
                raise ValueError ( "The copies of residue {} differ in atoms/protonation.".format ( name ) )
            r.atoms = [ ( a[0], a[1], np.asarray ( a[2], float ) ) for a in atoms ]
            r.library_name, r.kind, r.missing_heavy = name, "ligand", [ ]
            r.n_hydrogens = sum ( 1 for a in r.atoms if a[1] == 1 )
        analysis.ligands[name] = lig
    analysis.unsupported = [ ]
    if library_root not in analysis.extra_library_paths:
        analysis.extra_library_paths.insert ( 0, library_root )
    return dict ( analysis.ligands )
#=============================================================================
#  Typing
#=============================================================================
def parameter_set_path ( parameter_set = DEFAULT_PARAMETER_SET ):
    """ A parameter set name (pDynamo's opls/<name>) or a folder path. """
    if os.path.isdir ( parameter_set ):
        return os.path.abspath ( parameter_set )
    return os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "forceFields", "opls", parameter_set )

def type_atoms ( system, parameter_set = DEFAULT_PARAMETER_SET ):
    """ ( types [label or None per atom], charges, untyped [Atom] ) without
        building an MM model. """
    from pMolecule.MMModel.MMAtomTyper import MMAtomTyper
    from pMolecule.MMModel.MMModelError import MMModelError
    typer = MMAtomTyper.FromPath ( parameter_set_path ( parameter_set ) )
    types   = [ None ] * len ( system.atoms )
    charges = [ 0.0  ] * len ( system.atoms )
    untyped = set ( system.connectivity.atoms )
    typer.TypeBySequence ( system.sequence, types, charges, untyped )
    typer.TypeByPattern  ( system.connectivity, types, charges, untyped )
    return types, charges, sorted ( untyped, key = lambda a: a.index )

def typing_report ( system, parameter_set = DEFAULT_PARAMETER_SET ):
    """ { "ok", "n_atoms", "untyped" [paths], "untyped_by_residue" {residue: [atom names]},
          "types" Counter, "charge" (sum of typed charges), "formal_charge" } """
    types, charges, untyped = type_atoms ( system, parameter_set )
    by_residue = OrderedDict ( )
    for atom in untyped:
        by_residue.setdefault ( atom.parent.path, [ ] ).append ( atom.label )
    return { "ok"                 : len ( untyped ) == 0,
             "n_atoms"            : len ( types ),
             "untyped"            : [ a.path for a in untyped ],
             "untyped_by_residue" : by_residue,
             "types"              : Counter ( t for t in types if t is not None ),
             "charge"             : float ( sum ( charges ) ),
             "formal_charge"      : int ( sum ( getattr ( a, "formalCharge", 0 ) or 0 for a in system.atoms ) ) }


#=============================================================================
#  End to end
#=============================================================================
def copy_symmetry ( source, target ):
    """ The periodic cell (symmetry + parameters) of `source` on `target`
        (nothing when `source` has none). Returns True when a cell was set. """
    symmetry = getattr ( source, "symmetry", None )
    parameters = getattr ( source, "symmetryParameters", None )
    if symmetry is None or parameters is None: return False
    from pScientific.Symmetry import PeriodicBoundaryConditions
    crystal = symmetry.crystalSystem
    unique = crystal.GetUniqueSymmetryParameters ( parameters )
    target.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( crystal.__class__ ( ) )
    target.symmetryParameters = target.symmetry.MakeSymmetryParameters ( **unique )
    return True


def prepare_opls_system ( system, parameter_set = DEFAULT_PARAMETER_SET, work_dir = None, label = None,
                          analysis = None, define_nb_model = True, solvation = None, progress = None,
                          parametrize_unsupported = True, ligand_ph = 7.4, cm5_scale = None, xtb = None,
                          trust_connectivity = False, complete_truncated = False ):
    """ Loaded system -> prepared, typed system with an OPLS MM model.

        analysis: an edited analyze_system() result (protonation choices);
        solvation: None, or the options of solvate.solvate_coordinates()
        (shape, padding, box_size, neutralize, concentration, cation, anion, ...);
        progress: optional callable ( text ) for status messages;
        work_dir: where the per-system parameter set is written (default:
        a new temporary folder). Returns ( system, report ) where report =
        build report + "typing" + "parameters" (materialize report). The
        MM model is NOT defined when atoms stay untyped or parameters are
        missing; report["ok"] tells. """
    from pdynamo.opls import parameters
    from pMolecule.MMModel import MMModelOPLS
    from pMolecule.NBModel import NBModelCutOff
    say = progress or ( lambda text: None )
    if work_dir is None:
        work_dir = tempfile.mkdtemp ( prefix = "opls_set_" )
    os.makedirs ( work_dir, exist_ok = True )
    if analysis is None:
        say ( "Analyzing residues..." )
        analysis = analyze_system ( system )
    if analysis.unsupported and parametrize_unsupported:
        parametrize_ligands ( analysis, work_dir, ph = ligand_ph, cm5_scale = cm5_scale, progress = say, xtb = xtb,
                              trust_connectivity = trust_connectivity )
    say ( "Building hydrogens, termini and links..." )
    new_system, report = build_prepared_system ( analysis, label = label, complete_truncated = complete_truncated )
    # . [EN] 2026-10-17 BUG FIX (user's report: a system solvated in the Builder lost
    #   its box when OPLS was assigned): the rebuilt system keeps the input's periodic
    #   cell unless this call solvates it again (which sets its own cell)
    if solvation is None:
        copy_symmetry ( system, new_system )
    # . ligand total charges set by the user (not carried by the atoms' formal charges)
    correction = 0
    for r in analysis.residues:
        lig = analysis.ligands.get ( r.library_name ) if r.kind == "ligand" else None
        if lig is not None and len ( lig.atoms ) > 1:
            correction += lig.total_charge - lig.formal_charge_sum
    report["total_charge"] = report["formal_charge"] + correction
    if solvation is not None and not report["undefined_atoms"]:
        from pdynamo.opls import solvate
        say ( "Solvating..." )
        new_system, solvation_report = solvate.solvate_prepared_system ( new_system, label = label, ligands = dict ( analysis.ligands ),
                                                                         extra_library_paths = analysis.extra_library_paths,
                                                                         charge_correction = correction, **solvation )
        report["solvation"] = solvation_report
        report["atoms"] = solvation_report["atoms"]
        report["formal_charge"] = solvation_report["formal_charge"]
        report["total_charge"]  = solvation_report["formal_charge"] + correction
    if getattr ( analysis, "ligand_variants", None ):
        report["ligand_variants"] = dict ( analysis.ligand_variants )
    if analysis.ligands:
        from pdynamo.opls import parameters as _parameters
        parameter_set = _parameters.extend_parameter_set ( os.path.join ( work_dir, "base_with_ligands" ),
                                                           base = parameter_set, ligands = analysis.ligands.values ( ) )
        report["ligands"] = { name: lig.summary ( ) for name, lig in analysis.ligands.items ( ) }
    say ( "Typing atoms..." )
    typing = typing_report ( new_system, parameter_set )
    report["typing"] = typing
    report["ok"] = typing["ok"] and not report["undefined_atoms"]
    if not report["ok"]:
        return new_system, report
    types, charges, _ = type_atoms ( new_system, parameter_set )
    say ( "Writing the parameter set..." )
    database = parameters.build_class_database ( parameter_set )
    fallback = None
    if analysis.ligands:
        from pdynamo.opls import ligand as L
        fallback, notes = L.estimate_missing_terms ( new_system, types, database,
                                                     parameters.PDynamoSet ( parameters.pdynamo_set_path ( parameter_set ) ) )
        report["estimated_terms"] = notes
    set_dir = os.path.join ( work_dir, "parameters" )
    materialized = parameters.materialize_parameter_set ( new_system, types, set_dir, base = parameter_set,
                                                          database = database, fallback_terms = fallback )
    report["parameters"] = materialized
    if any ( materialized["missing"].values ( ) ):
        report["ok"] = False
        return new_system, report
    if analysis.ligands:
        from pdynamo.opls import quality
        try:
            report["quality"] = quality.assess ( new_system, types, charges, materialized, set ( analysis.ligands ) )
            import json
            with open ( os.path.join ( work_dir, "quality.json" ), "w" ) as handle:
                json.dump ( report["quality"], handle, indent = 1, default = str )
        except Exception as error:
            report["quality_error"] = str ( error )
    say ( "Building the MM model..." )
    new_system.DefineMMModel ( MMModelOPLS.WithParameterSet ( os.path.abspath ( set_dir ) ), log = None )
    if define_nb_model:
        new_system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
    report["parameter_folder"] = os.path.abspath ( set_dir )
    return new_system, report


#=============================================================================
#  Report
#=============================================================================
def report_markdown ( report, analysis = None, title = "OPLS system preparation" ):
    """ Human-readable report (Markdown) of prepare_opls_system(). """
    import datetime
    lines = [ "# {}".format ( title ), "", "Generated by EasyHybrid (pdynamo/opls), {}.".format (
              datetime.datetime.now ( ).strftime ( "%Y-%m-%d %H:%M" ) ), "" ]
    lines += [ "## Result", "",
               "- Status: **{}**".format ( "ready (OPLS MM model defined)" if report.get ( "ok" ) else "FAILED" ),
               "- Atoms: {} (input {}, hydrogens built {})".format ( report.get ( "atoms" ), report.get ( "input_atoms" ), report.get ( "added_hydrogens" ) ),
               "- Total charge: {:+d}{}".format ( report.get ( "total_charge", report.get ( "formal_charge", 0 ) ),
                    "" if report.get ( "total_charge", report.get ( "formal_charge" ) ) == report.get ( "formal_charge" )
                    else " (formal charges of the atoms: {:+d}; ligand charges set by the user)".format ( report.get ( "formal_charge", 0 ) ) ),
               "- Segments: {}".format ( ", ".join ( report.get ( "segments", [ ] ) ) ),
               "- Disulfide bridges: {}; chain breaks: {}".format ( report.get ( "disulfides", 0 ), report.get ( "chain_breaks", 0 ) ) ]
    if report.get ( "parameter_folder" ):
        lines.append ( "- Parameter set: `{}`".format ( report["parameter_folder"] ) )
    if analysis is not None:
        summary = analysis.summary ( )
        changed = { k: v for k, v in summary["variants"].items ( ) if v[1] != "library default" or v[0] }
        lines += [ "", "## Protonation states", "", "Library defaults (pH 7: ASP/GLU deprotonated, LYS/ARG protonated, HIS delta-protonated) except:", "" ]
        lines += [ "- {}: {} ({})".format ( k, v[0] or "library form", v[1] ) for k, v in changed.items ( ) ] or [ "- none" ]
        if summary["added_atoms"]:
            lines += [ "", "Atoms placed by geometry: " + ", ".join ( "{} {}".format ( k, "/".join ( v ) ) for k, v in summary["added_atoms"].items ( ) ) ]
        truncated = [ r.label for r in analysis.residues if r.variant_source == "truncated" ]
        if truncated:
            lines += [ "", "Truncated residues (incomplete side chains): " + ", ".join ( truncated ) ]
        rebuilt = [ "{} ({})".format ( r.label, " ".join ( r.rebuilt_atoms ) ) for r in analysis.residues if getattr ( r, "rebuilt_atoms", None ) ]
        if rebuilt:
            lines += [ "", "Rebuilt residues (missing heavy atoms placed geometrically -- minimize and check them): " + ", ".join ( rebuilt ) ]
    solvation = report.get ( "solvation" )
    if solvation:
        lines += [ "", "## Solvation", "",
                   "- Shape: {}".format ( solvation["shape"] ),
                   "- Cell: {}".format ( "{:.2f} x {:.2f} x {:.2f} A".format ( *solvation["cell"][:3] ) if solvation.get ( "cell" ) else "none (droplet, radius {:.2f} A)".format ( solvation["sphere"][1] ) ),
                   "- Waters (TIP3P): {}; removed by the solute: {}".format ( solvation["waters"], solvation["removed_by_solute"] ),
                   "- Ions: {} (solute charge {:+d})".format ( ", ".join ( "{} {}".format ( n, k ) for k, n in solvation["ions"].items ( ) ) or "none", solvation["solute_charge"] ) ]
    ligands = report.get ( "ligands" )
    if ligands:
        lines += [ "", "## Ligands (simplified parametrization)", "" ]
        for new, old in ( report.get ( "ligand_variants" ) or { } ).items ( ):
            lines.append ( "- {}: copy of {} with different atoms/bonds (edited), parametrized separately under this name.".format ( new, old ) )
        if report.get ( "ligand_variants" ): lines.append ( "" )
        for name, data in ligands.items ( ):
            lines.append ( "- **{}**: {} atoms, charge {:+d}, multiplicity {}, charges: {}, classes: {}".format (
                    name, data["atoms"], data["total_charge"], data.get ( "multiplicity", 1 ), data["charge_method"],
                    ", ".join ( "{} {}".format ( c, n ) for c, n in data["classes"].items ( ) ) ) )
            for note in data["notes"] + data["warnings"]:
                lines.append ( "  - {}".format ( note ) )
        estimated = report.get ( "estimated_terms", [ ] )
        lines += [ "", "### Estimated bonded terms ({})".format ( len ( estimated ) ), "" ]
        lines += [ "- " + e for e in estimated ] or [ "None: every term came from the OPLS database." ]
    parameters = report.get ( "parameters" )
    if parameters:
        def friendly ( source ):
            if source == "base":      return "pDynamo set, exact labels"
            if source.startswith ( "pdynamo:" ): return "pDynamo {} set, by class".format ( source.split ( ":", 1 )[1] )
            return { "oplsaal": "OPLS-AA/L (oplsaal.prm)", "oplsaa08": "OPLS-AA 2008 (oplsaa08.prm)",
                     "estimated": "estimated (see above)", "zero": "zero (no OPLS improper)" }.get ( source, source )
        lines += [ "", "## Parameter sources", "", "| Term | Sources |", "|---|---|" ]
        for term, counts in parameters["sources"].items ( ):
            if counts: lines.append ( "| {} | {} |".format ( term, ", ".join ( "{} {}".format ( friendly ( k ), v ) for k, v in counts.items ( ) ) ) )
        if parameters.get ( "zero_outofplane" ):
            lines.append ( "" )
            lines.append ( "Out-of-plane terms without parameters (K = 0): {}".format ( len ( parameters["zero_outofplane"] ) ) )
        missing = { t: v for t, v in parameters["missing"].items ( ) if v }
        if missing:
            lines += [ "", "**Missing parameters:**", "" ] + [ "- {}: {}".format ( t, ", ".join ( "-".join ( k ) for k in v ) ) for t, v in missing.items ( ) ]
    if report.get ( "quality" ):
        from pdynamo.opls import quality
        lines += quality.markdown ( report["quality"] )
    typing = report.get ( "typing" )
    if typing and not typing["ok"]:
        lines += [ "", "## Untyped atoms", "" ] + [ "- {}: {}".format ( r, ", ".join ( a ) ) for r, a in typing["untyped_by_residue"].items ( ) ]
    lines += [ "", "## Before production", "",
               "Minimize, then equilibrate (NVT, then NPT for periodic boxes). Ligand parameters are a simplified",
               "OPLS-AA/CM5 model: check the estimated terms above." ]
    return "\n".join ( lines ) + "\n"
