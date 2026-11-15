#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Test set of well-studied metalloproteins (RCSB PDB). Downloads the PDB files to
#  this folder, checks that the expected metals are there, and writes manifest.json
#  with the title, method, resolution and the metal sites read from the LINK records
#  (or, when there are none, from distances: heavy atoms within 2.8 A of the metal).
#
#      python3 tests/data/metalloproteins/fetch_metalloproteins.py
#

import json
import math
import os
import time
import urllib.request

HERE = os.path.dirname ( os.path.abspath ( __file__ ) )

# . id, protein, expected metal residue names, site in the literature, why it is in the set
SET = [
  ( "1ZNF", "Zinc finger Xfin-31 (NMR)",                    [ "ZN" ],        "Zn Cys2His2 (structural)", "classic structural Zn site; NMR model, caps ACE/NH2" ),
  ( "3KS3", "Carbonic anhydrase II, human (0.9 A)",         [ "ZN" ],        "Zn His3 + OH-/H2O (catalytic)", "the reference catalytic Zn site" ),
  ( "8TLN", "Thermolysin",                                   [ "ZN", "CA" ],  "Zn His2Glu + H2O; four Ca sites", "Zn protease with structural Ca" ),
  ( "4INS", "Insulin hexamer (2Zn)",                         [ "ZN" ],        "Zn His3 x2 (on the 3-fold axis)", "Zn at symmetry axis, multi-chain" ),
  ( "2C9V", "Cu,Zn superoxide dismutase 1, human",           [ "CU", "ZN" ],  "Cu His4 (bridging His); Zn His3Asp", "binuclear Cu/Zn with bridging imidazolate" ),
  ( "4AZU", "Azurin, Pseudomonas aeruginosa",                [ "CU" ],        "Cu(II) type 1: His2 Cys Met (+Gly O)", "blue copper, open-shell Cu(II)" ),
  ( "1PLC", "Plastocyanin, poplar (1.33 A)",                 [ "CU" ],        "Cu type 1: His2 Cys Met", "blue copper reference" ),
  ( "1IRO", "Rubredoxin, Clostridium pasteurianum (1.1 A)",  [ "FE" ],        "Fe Cys4 (high-spin Fe(II)/Fe(III))", "simplest Fe-S site, open shell" ),
  ( "1FDN", "Ferredoxin 2[4Fe-4S], C. acidurici",           [ "SF4" ],       "two [4Fe-4S] clusters, Cys4 each", "iron-sulfur clusters (cluster residues)" ),
  ( "1A6M", "Myoglobin, sperm whale, oxy (1.0 A)",           [ "HEM" ],       "heme Fe + proximal His + O2", "heme cofactor (porphyrin as one residue)" ),
  ( "1HRC", "Cytochrome c, horse heart",                     [ "HEC" ],       "heme c, Fe His/Met, thioether links to Cys", "covalently attached heme" ),
  ( "1CLL", "Calmodulin, human (Ca4)",                       [ "CA" ],        "four EF-hand Ca sites (Asp/Glu carboxylates, O)", "hard ion, O-only coordination" ),
  ( "5P21", "H-Ras p21 with GppNHp and Mg",                  [ "MG" ],        "Mg Ser/Thr + phosphate oxygens + waters", "Mg with nucleotide" ),
  ( "1UBQ", "Ubiquitin, human (no metal)",                  [ ],             "negative control: no metal", "checks that nothing is detected" ),
]

METALS = { "ZN", "CU", "FE", "MN", "CO", "NI", "MG", "CA", "NA", "K", "CD", "HG", "MO", "W" }


def fetch ( pdb_id ):
    path = os.path.join ( HERE, "{}.pdb".format ( pdb_id ) )
    if not os.path.exists ( path ):
        url = "https://files.rcsb.org/download/{}.pdb".format ( pdb_id )
        with urllib.request.urlopen ( url, timeout = 60 ) as response, open ( path, "wb" ) as handle:
            handle.write ( response.read ( ) )
        time.sleep ( 0.3 )
    return path


def parse ( path ):
    title, method, resolution, models = [ ], None, None, 0
    atoms, links = [ ], [ ]
    first_model = True
    for line in open ( path, errors = "replace" ):
        rec = line[:6].strip ( )
        if rec == "TITLE": title.append ( line[10:80].strip ( ) )
        elif rec == "EXPDTA": method = line[10:80].strip ( )
        elif rec == "REMARK" and line[7:10] == "  2" and "RESOLUTION." in line:
            try: resolution = float ( line.split ( "RESOLUTION." )[1].split ( )[0] )
            except ( ValueError, IndexError ): pass
        elif rec == "MODEL": models += 1
        elif rec == "ENDMDL": first_model = False
        elif rec == "LINK":
            links.append ( ( line[12:16].strip ( ), line[17:20].strip ( ), line[21], line[22:26].strip ( ),
                             line[42:46].strip ( ), line[47:50].strip ( ), line[51], line[52:56].strip ( ), line[73:78].strip ( ) ) )
        elif rec in ( "ATOM", "HETATM" ) and first_model:
            atoms.append ( { "name": line[12:16].strip ( ), "res": line[17:20].strip ( ), "chain": line[21], "seq": line[22:26].strip ( ),
                             "alt": line[16], "element": line[76:78].strip ( ).upper ( ) or line[12:14].strip ( ).upper ( ),
                             "xyz": ( float ( line[30:38] ), float ( line[38:46] ), float ( line[46:54] ) ), "het": rec == "HETATM" } )
    return " ".join ( title ), method, resolution, models, atoms, links


def sites ( atoms, links ):
    """ Metal sites: from LINK records involving a metal; else heavy atoms within 2.8 A. """
    out = [ ]
    metal_atoms = [ a for a in atoms if a["element"] in METALS and a["alt"] in ( " ", "A" ) ]
    for m in metal_atoms:
        key = ( m["name"], m["res"], m["chain"], m["seq"] )
        ligands = [ ]
        for l in links:
            for me, other in ( ( l[0:4], l[4:8] ), ( l[4:8], l[0:4] ) ):
                if ( me[0], me[1], me[2], me[3] ) == key:
                    ligands.append ( "{} {}{} {}".format ( other[1], other[2], other[3], other[0] ) )
        source = "LINK"
        if not ligands:
            source = "distance"
            for a in atoms:
                if a is m or a["element"] == "H" or a["alt"] not in ( " ", "A" ): continue
                d = math.dist ( a["xyz"], m["xyz"] )
                if d < 2.8: ligands.append ( "{} {}{} {} ({:.2f})".format ( a["res"], a["chain"], a["seq"], a["name"], d ) )
        out.append ( { "metal": "{} {}{}".format ( m["res"], m["chain"], m["seq"] ), "element": m["element"], "ligands": sorted ( set ( ligands ) ),
                       "source": source } )
    return out


def main ( ):
    manifest = [ ]
    for pdb_id, protein, expected, site, why in SET:
        entry = { "id": pdb_id, "protein": protein, "expected": expected, "site_literature": site, "why": why }
        try:
            path = fetch ( pdb_id )
            title, method, resolution, models, atoms, links = parse ( path )
            residues = { a["res"] for a in atoms }
            entry.update ( title = title, method = method, resolution = resolution, models = models or 1,
                           atoms_first_model = len ( atoms ), metal_sites = sites ( atoms, links ),
                           expected_found = all ( e in residues for e in expected ), file = os.path.basename ( path ) )
        except Exception as error:
            entry["error"] = str ( error )[:200]
        manifest.append ( entry )
        print ( "{:5} {:46} {:>6} {:>5}  sites {:2d}  expected found: {}".format ( pdb_id, protein[:46], ( entry.get ( "method" ) or "?" )[:6],
                entry.get ( "resolution" ) or "-", len ( entry.get ( "metal_sites", [ ] ) ), entry.get ( "expected_found", entry.get ( "error" ) ) ) )
    with open ( os.path.join ( HERE, "manifest.json" ), "w" ) as handle:
        json.dump ( { "created": time.strftime ( "%Y-%m-%d" ), "source": "RCSB PDB (https://www.rcsb.org)", "entries": manifest }, handle, indent = 1 )


if __name__ == "__main__":
    main ( )
