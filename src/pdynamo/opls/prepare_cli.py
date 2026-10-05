#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: OPLS system preparation from the command line.
#
#  Usage (from EasyHybrid3_dev/src):
#      python3 -m pdynamo.opls.prepare_cli  input.pdb  output.pkl  [--truncate-incomplete]
#              [--solvate cube|box|sphere] [--padding 10] [--salt 0.15] [--no-neutralize]
#              [--ph 7.4] [--cm5-scale 1.00] [--xtb /path/to/xtb]
#              [--ligand-charge RES=Q ...] [--ligand-mult RES=M ...]
#
#  Writes output.pkl (pDynamo system with the OPLS MM model and a cut-off NB
#  model, open it in EasyHybrid with Import System > pdynamo files), the
#  per-system parameter set, ligand files and prepare_opls_report.md in output_opls/
#

import os
import sys
import argparse


def main ( argv = None ):
    parser = argparse.ArgumentParser ( description = "Prepare a PDB system with OPLS (EasyHybrid)." )
    parser.add_argument ( "input" )
    parser.add_argument ( "output", help = "output .pkl" )
    parser.add_argument ( "--truncate-incomplete", action = "store_true",
                          help = "turn residues with missing heavy atoms into ALA/GLY" )
    parser.add_argument ( "--parameter-set", default = "protein" )
    parser.add_argument ( "--solvate", choices = ( "cube", "box", "sphere" ), default = None )
    parser.add_argument ( "--padding", type = float, default = 10.0 )
    parser.add_argument ( "--salt", type = float, default = 0.15, help = "salt concentration (mol/L)" )
    parser.add_argument ( "--no-neutralize", action = "store_true" )
    parser.add_argument ( "--ph", type = float, default = 7.4, help = "pH for ligand hydrogens" )
    parser.add_argument ( "--cm5-scale", type = float, default = None )
    parser.add_argument ( "--xtb", default = None )
    parser.add_argument ( "--ligand-charge", action = "append", default = [ ], metavar = "RES=Q",
                          help = "total charge of ligand RES (QC charge calculation)" )
    parser.add_argument ( "--ligand-mult", action = "append", default = [ ], metavar = "RES=M",
                          help = "spin multiplicity of ligand RES" )
    args = parser.parse_args ( argv )

    from pBabel import ImportSystem
    from pCore import Pickle
    from pdynamo.opls import prep

    system   = ImportSystem ( args.input, log = None )
    analysis = prep.analyze_system ( system )
    summary  = analysis.summary ( )
    print ( "Residues      :", summary["kinds"] )
    print ( "Segments      :", summary["segments"] )
    print ( "Chain breaks  :", summary["chain_breaks"] or "none" )
    print ( "Disulfides    :", summary["disulfides"] or "none" )
    print ( "Added atoms   :", summary["added_atoms"] or "none" )
    changed = { k: v for k, v in summary["variants"].items ( ) if v[1] != "library default" }
    print ( "Protonation   : library defaults (pH 7){}".format ( ", except " + str ( changed ) if changed else "" ) )
    if summary["unsupported"]:
        for name, entry in prep.ligand_preview ( analysis, ph = args.ph ).items ( ):
            print ( "Ligand        : {} x{} ({} heavy atoms), perceived charge {}".format (
                    name, entry["copies"], entry["heavy"], "?" if entry["charge"] is None else "{:+d}".format ( entry["charge"] ) ) )
    for option, key in ( ( args.ligand_charge, "charge" ), ( args.ligand_mult, "multiplicity" ) ):
        for item in option:
            name, _, value = item.partition ( "=" )
            analysis.ligand_options.setdefault ( name.strip ( ), { } )[key] = int ( value )
    if summary["missing_heavy"]:
        if not args.truncate_incomplete:
            print ( "\nMissing heavy atoms:", summary["missing_heavy"], "\n(use --truncate-incomplete)" )
            return 1
        for i, r in enumerate ( analysis.residues ):
            if r.missing_heavy:
                print ( "Truncated     :", r.label, "->", prep.truncate_residue ( analysis, i ) )

    base     = os.path.splitext ( args.output )[0]
    work_dir = base + "_opls"
    label    = os.path.basename ( base )
    solvation = None
    if args.solvate:
        solvation = { "shape": args.solvate, "padding": args.padding, "concentration": args.salt,
                      "neutralize": not args.no_neutralize }
    new_system, report = prep.prepare_opls_system ( system, parameter_set = args.parameter_set,
                                                     work_dir = work_dir, label = label, analysis = analysis,
                                                     solvation = solvation, ligand_ph = args.ph,
                                                     cm5_scale = args.cm5_scale, xtb = args.xtb,
                                                     progress = lambda text: print ( "  ..", text ) )
    print ( "\nAtoms         : {} (input {}, hydrogens built {})".format ( report["atoms"], report["input_atoms"], report["added_hydrogens"] ) )
    print ( "Total charge  : {:+d}".format ( report.get ( "total_charge", report["formal_charge"] ) ) )
    if not report["ok"]:
        print ( "FAILED -- untyped:", dict ( report["typing"]["untyped_by_residue"] ) )
        if "parameters" in report:
            print ( "missing parameters:", { t: v for t, v in report["parameters"]["missing"].items ( ) if v } )
        return 1
    sources = { t: dict ( c ) for t, c in report["parameters"]["sources"].items ( ) if c }
    print ( "Parameters    :", work_dir, sources )
    print ( "MM charge     : {:.4f}".format ( report["typing"]["charge"] ) )
    print ( "Energy        : {:.3f} kcal/mol".format ( new_system.Energy ( log = None ) ) )
    Pickle ( args.output, new_system )
    with open ( os.path.join ( work_dir, "prepare_opls_report.md" ), "w" ) as handle:
        handle.write ( prep.report_markdown ( report, analysis ) )
    print ( "Written       :", args.output, "and", os.path.join ( work_dir, "prepare_opls_report.md" ) )
    return 0


if __name__ == "__main__":
    sys.exit ( main ( ) )
