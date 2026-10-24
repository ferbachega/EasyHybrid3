#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: DYFF residue library -- precomputed charges, bonds/angles and torsions
#          of amino acids applied to a protein
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      The library (src/pdynamo/data/dyff_residue_library.json, made by
#      fraglib/fraglib/residue_library.py from capped ACE-X-NME models with
#      xTB) holds, per residue template: its atoms (element, bonded
#      neighbours, backbone role), CM5 charges (integer residue total),
#      Seminario bonds/angles and fitted torsions, the junction terms written
#      with references to the neighbours ("p:C", "p:O", "p:CA", "n:N", "n:H",
#      "n:CA").
#
#      Matching: every residue of the system is compared with the templates
#      of its family (residue name, histidine/protonation aliases) by GRAPH
#      ISOMORPHISM (elements + bonds, networkx): the protonation variant is
#      recognised from the atoms actually present (HID/HIE/HIP, ASP/ASH, ...)
#      and atom names do not matter. Residues that match no template
#      (termini, disulfide-bonded Cys, ligands, metals, ...) are returned as
#      "unmatched" atoms -- also residues with a bond leaving the residue from
#      an atom other than the backbone N and C (disulfide Cys, a side chain
#      bonded to a metal or a ligand) -- their charges come from elsewhere (CM5 on a
#      capped subsystem, see dyff_parameters.residue_library_charges_for_rest).
#
#      Output: DYFF per-term edits in system indices ("bond_terms",
#      "angle_terms", "dihedral_bonds"; a term written by two neighbouring
#      residues -- the peptide junction -- is averaged) and a charge per
#      matched atom.
#

import json
import os
from collections import defaultdict

#  2026-10-08: the default is the force-matched library (fraglib/run_force_matching.py:
#  angle scales and torsion series of the internal residues fitted to DFT forces and
#  energies of the SPICE dipeptides, training half; on the test half heavy-atom force
#  RMSE 132.7 -> 91.1 kJ/mol/A, OPLS-AA 96.7). It is built on the calibrated library
#  (Seminario constants calibrated to DFT by environment); termini and CYX, absent from
#  SPICE, keep the calibrated values. The xTB-only and calibrated files are kept for
#  comparison (fraglib/reports/fraglib_status.md sections 13-14).
DEFAULT_PATH = os.path.join ( os.path.dirname ( os.path.abspath ( __file__ ) ), "data", "dyff_residue_library_fm.json" )
ALIASES = { "HIS": "HIS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
            "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "CYM": "CYS", "CYX": "CYS" }


def load_library ( path = None ):
    with open ( path or DEFAULT_PATH ) as handle:
        library = json.load ( handle )
    library["templates"] = { k: v for k, v in library["templates"].items ( ) if "error" not in v }
    return library


def _template_graph ( template ):
    import networkx as nx
    g = nx.Graph ( )
    for k, atom in enumerate ( template["atoms"] ):
        g.add_node ( k, element = atom["element"] )
    for k, atom in enumerate ( template["atoms"] ):
        for n in atom["neighbors"]: g.add_edge ( k, n )
    return g


def match_residues ( residues, symbols, neighbours, library ):
    """ residues: [ ( residue label, generic name, [ atom indices ] ) ];
        symbols: element of every system atom; neighbours: [ set of bonded atoms ].
        Returns ( matches, unmatched ): matches = [ { "residue", "template",
        "map" ( template atom -> system atom ), "junction" ( ref -> system atom ) } ]. """
    import networkx as nx
    from networkx.algorithms import isomorphism
    graphs = { name: _template_graph ( t ) for name, t in library["templates"].items ( ) }
    by_family = defaultdict ( list )
    for name, t in library["templates"].items ( ): by_family[t["family"]].append ( name )
    matches, unmatched = [ ], [ ]
    for label, generic, atoms in residues:
        family = ALIASES.get ( generic.upper ( ), generic.upper ( ) )
        members = set ( atoms )
        g = nx.Graph ( )
        for a in atoms: g.add_node ( a, element = symbols[a] )
        for a in atoms:
            for b in neighbours[a]:
                if b in members: g.add_edge ( a, b )
        found = None
        for name in by_family.get ( family, [ ] ):
            t = graphs[name]
            if t.number_of_nodes ( ) != g.number_of_nodes ( ) or t.number_of_edges ( ) != g.number_of_edges ( ): continue
            gm = isomorphism.GraphMatcher ( t, g, node_match = lambda x, y: x["element"] == y["element"] )
            if not gm.is_isomorphic ( ): continue
            mapping = dict ( gm.mapping )
            # . the internal graph is not enough: only the backbone N and C may be bonded
            #   outside the residue (a disulfide-bonded Cys has the internal graph of a
            #   thiolate CYM -- its S is bonded to another residue)
            #   ("X": the atom of a CYX template bonded to the partner cysteine's S, x:S)
            backbone = { mapping[k] for k, atom in enumerate ( library["templates"][name]["atoms"] ) if atom["role"] in ( "N", "C", "X" ) }
            if any ( b not in members for a in atoms if a not in backbone for b in neighbours[a] ): continue
            xs = [ mapping[k] for k, atom in enumerate ( library["templates"][name]["atoms"] ) if atom["role"] == "X" ]
            if any ( len ( [ b for b in neighbours[x] if b not in members and symbols[b] == "S" ] ) != 1 for x in xs ): continue
            found = ( name, mapping ); break
        if found is None:
            unmatched.append ( ( label, generic, atoms ) ); continue
        name, mapping = found
        template = library["templates"][name]
        role = { atom["role"]: mapping[k] for k, atom in enumerate ( template["atoms"] ) if atom["role"] }
        junction = { }
        N, C = role.get ( "N" ), role.get ( "C" )
        pc = [ b for b in neighbours[N] if b not in members and symbols[b] == "C" ] if N is not None else [ ]
        if len ( pc ) == 1:
            junction["p:C"] = pc[0]
            po = [ b for b in neighbours[pc[0]] if symbols[b] == "O" and len ( neighbours[b] ) == 1 ]
            pca = [ b for b in neighbours[pc[0]] if symbols[b] == "C" and b != N ]
            if po: junction["p:O"] = po[0]
            if pca: junction["p:CA"] = pca[0]
        nn = [ b for b in neighbours[C] if b not in members and symbols[b] == "N" ] if C is not None else [ ]
        if len ( nn ) == 1:
            junction["n:N"] = nn[0]
            nh = [ b for b in neighbours[nn[0]] if symbols[b] == "H" ]
            nca = [ b for b in neighbours[nn[0]] if symbols[b] == "C" and b != C ]
            if len ( nh ) == 1: junction["n:H"] = nh[0]
            if nca: junction["n:CA"] = min ( nca, key = lambda b: len ( neighbours[b] ) != 4 )
        for k, atom in enumerate ( template["atoms"] ):
            if atom["role"] == "X":
                junction["x:S"] = [ b for b in neighbours[mapping[k]] if b not in members and symbols[b] == "S" ][0]
        matches.append ( { "residue": label, "template": name, "map": mapping, "junction": junction } )
    return matches, unmatched


def library_edits ( matches, library ):
    """ ( per-term DYFF edits in system indices, { atom: charge } ). Terms whose
        neighbour references are missing (chain ends) are left out; a term given
        by two residues is averaged. """
    bonds, angles, torsions, charges = defaultdict ( list ), defaultdict ( list ), defaultdict ( list ), { }
    for m in matches:
        t = library["templates"][m["template"]]
        def atom ( ref ):
            kind, key = ref.split ( ":" )
            return m["map"].get ( int ( key ) ) if kind == "r" else m["junction"].get ( ref )
        for k, a in enumerate ( t["atoms"] ):
            charges[m["map"][k]] = float ( a["charge"] )
        for ra, rb, b0, k in t["bond_terms"]:
            a, b = atom ( ra ), atom ( rb )
            if a is None or b is None: continue
            bonds["{}-{}".format ( min ( a, b ), max ( a, b ) )].append ( ( b0, k ) )
        for ra, rb, rc, pairs in t["angle_terms"]:
            a, b, c = atom ( ra ), atom ( rb ), atom ( rc )
            if None in ( a, b, c ): continue
            angles["{}-{}-{}".format ( min ( a, c ), b, max ( a, c ) )].append ( pairs )
        for ra, rb, pairs in t["dihedral_bonds"]:
            a, b = atom ( ra ), atom ( rb )
            if a is None or b is None: continue
            torsions["{}-{}".format ( min ( a, b ), max ( a, b ) )].append ( pairs )
    def average_pairs ( lists ):
        total = defaultdict ( float )
        for pairs in lists:
            for c, n in pairs: total[int ( n )] += float ( c )
        return [ [ total[n] / len ( lists ), n ] for n in sorted ( total ) ]
    edits = { "bond_terms": { k: [ sum ( v[0] for v in vals ) / len ( vals ), sum ( v[1] for v in vals ) / len ( vals ) ]
                              for k, vals in bonds.items ( ) },
              "angle_terms": { k: average_pairs ( vals ) for k, vals in angles.items ( ) },
              "dihedral_bonds": { k: average_pairs ( vals ) for k, vals in torsions.items ( ) } }
    return edits, charges


def residues_of_system ( system ):
    """ [ ( path, generic name, [ atom indices ] ) ] from a pDynamo system's sequence. """
    out = [ ]
    if system.sequence is None: return out
    for entity in system.sequence.children:
        for component in entity.children:
            out.append ( ( component.path, component.genericLabel, [ a.index for a in component.children ] ) )
    return out


def system_graph ( system ):
    from pScientific import PeriodicTable
    symbols = [ PeriodicTable.Symbol ( a.atomicNumber ) for a in system.atoms ]
    neighbours = [ set ( ) for _ in system.atoms ]
    for i, j in system.connectivity.bondIndices:
        neighbours[i].add ( j ); neighbours[j].add ( i )
    return symbols, neighbours


#=============================================================================
#  Charges for the atoms the library does not cover
#=============================================================================
def xtb_cm5 ( symbols, xyz, charge, solvent = "water", threads = 4, multiplicity = 1, log_dir = None ):
    """ GFN1-xTB CM5 charges (implicit solvent by default; multiplicity > 1: --uhf).
        log_dir: the input geometry and the xtb output are kept there (param_logs). """
    import shutil, subprocess, tempfile
    import numpy as np
    from pdynamo.opls.ligand import find_xtb
    tmp = tempfile.mkdtemp ( prefix = "cm5_" )
    try:
        with open ( os.path.join ( tmp, "m.xyz" ), "w" ) as handle:
            handle.write ( "{}\n\n".format ( len ( symbols ) ) )
            for s, p in zip ( symbols, xyz ): handle.write ( "{} {:.6f} {:.6f} {:.6f}\n".format ( s, *p ) )
        command = [ find_xtb ( ), "m.xyz", "--gfn", "1", "--chrg", str ( int ( charge ) ) ] + ( [ "--alpb", solvent ] if solvent else [ ] )
        if multiplicity > 1: command += [ "--uhf", str ( int ( multiplicity ) - 1 ) ]
        r = subprocess.run ( command, cwd = tmp, capture_output = True, text = True, env = dict ( os.environ, OMP_NUM_THREADS = str ( threads ) ) )
        from pdynamo import param_logs
        param_logs.keep ( tmp, log_dir, [ "m.xyz" ], { "xtb_cm5.out": "$ " + " ".join ( command[1:] ) + "\n" + r.stdout + r.stderr } )
        start = r.stdout.find ( "Mulliken/CM5 charges" )
        if start < 0: raise RuntimeError ( "xtb gave no CM5 charges: " + ( r.stdout + r.stderr )[-400:] )
        return np.array ( [ float ( l.split ( )[2] ) for l in r.stdout[start:].splitlines ( )[1:len ( symbols ) + 1] ] )
    finally:
        shutil.rmtree ( tmp, ignore_errors = True )


def rest_charges ( atomic_numbers, coordinates, bonds, formal, atoms, shell = 2, solvent = "water", threads = 4, multiplicity = 1,
                   log_dir = None ):
    """ CM5 charges for `atoms` (the ones the library did not cover), computed on a
        capped subsystem around them (pdynamo/subsystem.py); every connected group of
        those atoms gets an integer total (the nearest integer of its CM5 sum).
        bonds: [ ( a, b, order, aromatic ) ] of the whole system. Returns { atom: charge }. """
    import numpy as np
    from pScientific import PeriodicTable
    from pdynamo import subsystem as SS
    atoms = sorted ( set ( atoms ) )
    if not atoms: return { }
    sub = SS.build_subsystem ( atomic_numbers, coordinates, bonds, atoms, shell = shell, formal_charges = formal )
    symbols = [ PeriodicTable.Symbol ( z ) for z in sub["atomic_numbers"] ]
    q = xtb_cm5 ( symbols, sub["coordinates"], sub["charge"], solvent = solvent, threads = threads, multiplicity = multiplicity,
                  log_dir = log_dir )
    local = { p: k for k, p in enumerate ( sub["parent_index"] ) if p >= 0 }
    charges = { a: float ( q[local[a]] ) for a in atoms }
    # . connected groups of the rest atoms -> integer totals
    members = set ( atoms )
    nb = defaultdict ( set )
    for a, b, *_ in bonds:
        if a in members and b in members: nb[a].add ( b ); nb[b].add ( a )
    seen = set ( )
    for a in atoms:
        if a in seen: continue
        group, stack = [ ], [ a ]; seen.add ( a )
        while stack:
            u = stack.pop ( ); group.append ( u )
            for v in nb[u]:
                if v not in seen: seen.add ( v ); stack.append ( v )
        total = sum ( charges[u] for u in group )
        shift = ( round ( total ) - total ) / len ( group )
        for u in group: charges[u] += shift
    return charges


#=============================================================================
#  Whole workflow for a pDynamo system (System > Add > Force Field > DYFF + residue library)
#=============================================================================
#  Water and monatomic ions (2026-10-08): a system solvated by "Prepare OPLS System"
#  would otherwise send every water molecule to the xTB treatment of the atoms outside
#  the library (thousands of identical jobs). They get fixed charges instead: TIP3P
#  for water (the charges of EasyHybrid's OPLS water) and the formal charge for ions;
#  the water oxygen also gets the TIP3P Lennard-Jones (dyff_edits.water_lj). This is a
#  flexible TIP3P-charge water on DYFF bonded terms -- NOT a validated water model.
WATER_CHARGES = { "O": -0.834, "H": 0.417 }
ION_CHARGES = { "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Mg": 2, "Ca": 2, "Zn": 2, "Mn": 2, "Fe": 2, "Cu": 2, "Ni": 2, "Co": 2,
                "F": -1, "Cl": -1, "Br": -1, "I": -1 }

def solvent_and_ions ( symbols, neighbours, atoms, formal = None ):
    """ ( water atoms, ion atoms, { atom: charge } ) among `atoms`: water = an O bonded
        to exactly two H and nothing else (with those H); ion = an atom with no bond. """
    atoms = set ( atoms )
    water, ions, charges = [ ], [ ], { }
    for a in sorted ( atoms ):
        if symbols[a] == "O" and len ( neighbours[a] ) == 2 and all ( symbols[h] == "H" and len ( neighbours[h] ) == 1 and h in atoms for h in neighbours[a] ):
            group = [ a ] + sorted ( neighbours[a] )
            water += group
            for i in group: charges[i] = WATER_CHARGES[symbols[i]]
        elif not neighbours[a] and symbols[a] != "H":
            q = int ( formal[a] ) if formal is not None and formal[a] else ION_CHARGES.get ( symbols[a] )
            if q is not None:
                ions.append ( a ); charges[a] = float ( q )
    return water, ions, charges


def plan ( system, library = None ):
    """ Matching only (fast): { matches, unmatched, edits, charges, rest_atoms, summary };
        water and monatomic ions get fixed charges and are not "rest" atoms. """
    from collections import Counter
    library = library or load_library ( )
    symbols, neighbours = system_graph ( system )
    matches, unmatched = match_residues ( residues_of_system ( system ), symbols, neighbours, library )
    edits, charges = library_edits ( matches, library )
    covered = set ( charges )
    rest = [ i for i in range ( len ( system.atoms ) ) if i not in covered ]
    formal = [ int ( getattr ( a, "formalCharge", 0 ) or 0 ) for a in system.atoms ]
    water, ions, fixed = solvent_and_ions ( symbols, neighbours, rest, formal )
    charges.update ( fixed )
    rest = [ i for i in rest if i not in fixed ]
    water_residues = { label for label, generic, atoms in unmatched if atoms and all ( a in fixed for a in atoms ) }
    unmatched = [ u for u in unmatched if u[0] not in water_residues ]
    summary = { "residues": len ( matches ) + len ( unmatched ), "matched": len ( matches ),
                "templates": dict ( Counter ( m["template"] for m in matches ) ),
                "unmatched": [ "{} {}".format ( g, l.split ( ":" )[-1] ) for l, g, a in unmatched ], "rest_atoms": len ( rest ),
                "waters": len ( water ) // 3, "ions": len ( ions ) }
    return { "matches": matches, "unmatched": unmatched, "edits": edits, "charges": charges, "rest_atoms": rest, "summary": summary }


def system_bonds ( system ):
    out = [ ]
    for b in system.connectivity.bonds:
        order = { "Single": 1, "Double": 2, "Triple": 3 }.get ( getattr ( b.type, "name", "Single" ), 1 )
        out.append ( ( b.node1.index, b.node2.index, order, bool ( getattr ( b, "isAromatic", False ) ) ) )
    return out


def compute_rest ( system, rest_atoms, solvent = "water", threads = 4 ):
    """ CM5 charges for the atoms outside the library (worker-thread safe). """
    import numpy as np
    xyz = np.array ( [ [ system.coordinates3[a, k] for k in range ( 3 ) ] for a in range ( len ( system.atoms ) ) ] )
    formal = [ int ( getattr ( a, "formalCharge", 0 ) or 0 ) for a in system.atoms ]
    return rest_charges ( [ a.atomicNumber for a in system.atoms ], xyz, system_bonds ( system ), formal, rest_atoms,
                          solvent = solvent, threads = threads )


# . ligand charge models of the DYFF dialogs (ligand_charges) -> pdynamo/esp_charges.py level
RESP_LEVELS = { "resp-xtb": "GFN2-xTB", "resp-hf": "HF/6-31G*" }
LIGAND_CHARGE_MODELS = ( ( "cm5", "CM5 (GFN1-xTB), default" ), ( "resp-xtb", "RESP on the GFN2-xTB potential (fast)" ),
                         ( "resp-hf", "RESP on the HF/6-31G* potential (AMBER level; ORCA when installed, else pDynamo up to ~22 atoms)" ) )


def resp_group_charges ( z, xyz, bonds, formal, atoms, level, multiplicity_of = None, log_dir = None, nprocs = None ):
    """ Two-stage RESP charges (pdynamo/esp_charges.py, gas phase) for every connected group
        of `atoms` that is a whole molecule: no covalent bond to atoms outside the group and
        no metal. Groups bonded to the rest of the system (termini, modified residues,
        covalent ligands) and metal sites keep their CM5 charges.
        log_dir: a folder resp_<first atom> per molecule with the original xtb/ORCA outputs;
        nprocs: ORCA cores.
        Returns ( { atom: charge }, [ { "atoms", "status", ... } ] ). """
    import numpy as np
    from pdynamo import esp_charges as EC
    from pdynamo.subsystem import _METALS
    multiplicity_of = multiplicity_of or ( lambda group: 1 )
    members = set ( atoms )
    nb = defaultdict ( set )
    for a, b, *_ in bonds: nb[a].add ( b ); nb[b].add ( a )
    groups, seen = [ ], set ( )
    for a in sorted ( members ):
        if a in seen: continue
        group, stack = [ ], [ a ]; seen.add ( a )
        while stack:
            u = stack.pop ( ); group.append ( u )
            for v in nb[u]:
                if v in members and v not in seen: seen.add ( v ); stack.append ( v )
        groups.append ( sorted ( group ) )
    out, info = { }, [ ]
    for group in groups:
        gset = set ( group )
        entry = { "atoms": len ( group ), "members": [ int ( u ) for u in group ] }
        if any ( v not in gset for u in group for v in nb[u] ):
            entry["status"] = "CM5 kept: bonded to the rest of the system"
        elif any ( z[u] in _METALS for u in group ):
            entry["status"] = "CM5 kept: metal"
        elif len ( group ) == 1:
            out[group[0]] = float ( formal[group[0]] ); entry["status"] = "ion: formal charge"
        else:
            local = { a: k for k, a in enumerate ( group ) }
            q = int ( sum ( formal[u] for u in group ) )
            try:
                r = EC.compute_esp_charges ( [ z[u] for u in group ], [ list ( xyz[u] ) for u in group ],
                                             [ ( local[a], local[b] ) for a, b, *_ in bonds if a in gset and b in gset ],
                                             q, multiplicity_of ( gset ), level,
                                             log_dir = os.path.join ( log_dir, "resp_{}".format ( group[0] + 1 ) ) if log_dir else None,
                                             nprocs = nprocs )
                for u, c in zip ( group, r["resp2"]["charges"] ): out[u] = float ( c )
                entry.update ( status = "RESP-2 ({})".format ( level ), charge = q, rrms = r["resp2"]["rrms"],
                               seconds = round ( r["seconds"], 1 ), engine = r.get ( "engine" ), level = level,
                               multiplicity = int ( multiplicity_of ( gset ) ), points = r.get ( "points" ),
                               dipole_qm = r.get ( "qc_dipole" ), dipole_resp = r["resp2"].get ( "dipole" ),
                               esp_charges = [ float ( c ) for c in r["esp"]["charges"] ],
                               log_dir = "resp_{}".format ( group[0] + 1 ) if log_dir else None )
            except Exception as error:
                entry["status"] = "CM5 kept: RESP failed ({})".format ( str ( error ).split ( "\n" )[0][:120] )
        info.append ( entry )
    return out, info


def rest_parameters ( system, rest_atoms, solvent = "water", threads = 4, hessian = True, torsions = True, shell = 2, say = None,
                      method = "gfn2", site_states = None, scan_torsions = False, known_charges = None, known_edits = None,
                      total_charge = None, multiplicity = None, scan_fragments = True, scan_gas = True,
                      ligand_charges = "cm5", log_dir = None ):
    """ Parameters for the atoms outside the library (ligands, cofactors, termini,
        modified residues), on ONE capped subsystem around them (pdynamo/subsystem.py):
          charges   GFN1-xTB CM5, integer total per connected group (as rest_charges);
          bonds/angles  modified Seminario on the GFN2-xTB Hessian of the subsystem;
          torsions  the DFT torsion library (pdynamo/bond_library.py), looked up with the
                    codes of the WHOLE system (the caps would change the environment).
        Only terms whose atoms are all rest atoms are returned (system indices); the
        torsion library was validated on top of Seminario bonds/angles, so torsions
        are taken only together with them. Worker-thread safe.
        method: "gfn2" (GFN2-xTB, implicit solvent allowed) or "gxtb" (g-xTB, gas phase:
        tblite has no implicit solvent) for the Hessians; the charges stay GFN1-xTB CM5.
        site_states: { metal atom index: { "multiplicity" } } (dyff_prep.site_groups):
        the spin of the groups containing those metals (unpaired electrons added up
        when a group / the charge calculation holds several open-shell metals).
        scan_torsions (2026-10-21): the flexible rotatable bonds the torsion library does
        not have (flexible_scan_targets: no rings, multiple or amide-like bonds, no
        CH3/NH2/NH3+ rotors) get a relaxed xTB scan and a fitted per-bond torsion
        (scan_flexible_torsions); known_charges / known_edits: the library's charges
        and terms of the other residues, used in those fits.
        log_dir (2026-10-23): the original xtb/ORCA inputs and outputs of every step are kept
        there (charges/, hessian/, torsions/), for the parametrization report
        (pdynamo/param_report.py); report["timings"] has the wall time of each step.
        Returns { "charges", "edits": { bond_terms, angle_terms, dihedral_bonds }, "report" }. """
    import numpy as np
    from pScientific import PeriodicTable
    from pdynamo import subsystem as SS
    import time
    say = say or ( lambda text: None )
    rest = sorted ( set ( rest_atoms ) )
    edits = { "bond_terms": { }, "angle_terms": { }, "dihedral_bonds": { } }
    report = { "rest_atoms": len ( rest ), "bonds": 0, "angles": 0, "torsions": 0, "rotatable": 0, "levels": { }, "skipped": 0,
               "timings": { } }
    sub_log = lambda name: os.path.join ( log_dir, name ) if log_dir else None
    clock = [ time.time ( ) ]
    def lap ( name ):
        now = time.time ( ); report["timings"][name] = round ( now - clock[0], 2 ); clock[0] = now
    if not rest: return { "charges": { }, "edits": edits, "report": report }
    n = len ( system.atoms )
    z = [ a.atomicNumber for a in system.atoms ]
    xyz = np.array ( [ [ system.coordinates3[a, k] for k in range ( 3 ) ] for a in range ( n ) ] )
    formal = [ int ( getattr ( a, "formalCharge", 0 ) or 0 ) for a in system.atoms ]
    bonds = system_bonds ( system )
    site_states = site_states or { }
    # . 2026-10-21, user's request: charge and multiplicity chosen in the dialog.
    #   total_charge: when it differs from the sum of the formal charges, the difference
    #   goes on the metal (one metal among the rest atoms) or on the first rest atom,
    #   for the xTB calculations only. multiplicity: the groups without a metal site
    #   (metal sites carry theirs in site_states, set by the caller).
    if total_charge is not None and int ( total_charge ) != sum ( formal ):
        offset = int ( total_charge ) - sum ( formal )
        rest_metals = [ a for a in rest if a in site_states ]
        target = rest_metals[0] if len ( rest_metals ) == 1 else rest[0]
        formal = list ( formal ); formal[target] += offset
        report["charge_offset"] = { "atom": int ( target ), "offset": int ( offset ) }
    def multiplicity_of ( atoms ):
        if multiplicity is not None and not any ( m in atoms for m in site_states ):
            return int ( multiplicity )
        unpaired = sum ( int ( st.get ( "multiplicity", 1 ) ) - 1 for m, st in site_states.items ( ) if m in atoms )
        return 1 + unpaired
    say ( "CM5 charges of the {} atoms outside the library...".format ( len ( rest ) ) )
    charges = rest_charges ( z, xyz, bonds, formal, rest, shell = shell, solvent = solvent, threads = threads,
                             multiplicity = multiplicity_of ( set ( rest ) ), log_dir = sub_log ( os.path.join ( "charges", "cm5" ) ) )
    report["cm5_charges"] = { int ( a ): float ( q ) for a, q in charges.items ( ) }
    report["charge_settings"] = { "cm5": "GFN1-xTB CM5" + ( ", ALPB({})".format ( solvent ) if solvent else ", gas phase" ),
                                  "ligand_charges": ligand_charges or "cm5" }
    lap ( "cm5" )
    # . 2026-10-23 (user decision): RESP charges for the ligands instead of CM5
    if ligand_charges and ligand_charges != "cm5":
        level = RESP_LEVELS[ligand_charges]
        say ( "RESP charges ({}) of the ligands...".format ( level ) )
        resp, info = resp_group_charges ( z, xyz, bonds, formal, rest, level, multiplicity_of,
                                          log_dir = sub_log ( "charges" ), nprocs = max ( 1, threads ) )
        charges.update ( resp )
        report["ligand_charges"] = info
        lap ( "resp" )
    report["method"] = method
    if hessian:
        import concurrent.futures
        from pdynamo import seminario as SM
        # . one capped subsystem per connected group of rest atoms (a terminus, a pair of
        #   bridged cysteines, a ligand): the numerical xTB Hessian costs ~N^2, one
        #   subsystem around ALL of them (crambin: 369 atoms) took tens of minutes
        members = set ( rest )
        nb = defaultdict ( set )
        for a, b, *_ in bonds:
            if a in members and b in members: nb[a].add ( b ); nb[b].add ( a )
        groups, seen = [ ], set ( )
        for a in rest:
            if a in seen: continue
            group, stack = [ ], [ a ]; seen.add ( a )
            while stack:
                u = stack.pop ( ); group.append ( u )
                for v in nb[u]:
                    if v not in seen: seen.add ( v ); stack.append ( v )
            groups.append ( sorted ( group ) )
        # . DFT/xTB calibration of the Seminario constants by environment (pdynamo/bond_library.py)
        calibration = None
        try:
            from pdynamo import bond_library as BL
            if BL.load_hessian_library ( ) is not None: calibration = BL.HessianCalibration ( z, bonds, formal, xyz )
        except Exception:
            calibration = None
        def one ( group ):
            sub = SS.build_subsystem ( z, xyz, bonds, group, shell = shell, formal_charges = formal )
            symbols = [ PeriodicTable.Symbol ( x ) for x in sub["atomic_numbers"] ]
            hdir = sub_log ( os.path.join ( "hessian", "group_{}".format ( group[0] + 1 ) ) )
            x_opt, h, _ = SM.xtb_hessian ( symbols, sub["coordinates"], charge = sub["charge"], threads = 1, solvent = solvent,
                                           multiplicity = multiplicity_of ( set ( group ) ), gfn = "gxtb" if method == "gxtb" else 2,
                                           log_dir = hdir )
            if hdir:
                from pdynamo import param_logs
                param_logs.keep ( None, hdir, texts = { "parent_atoms.txt": "subsystem atom -> system atom (1-based, 0 = cap)\n" +
                                  "\n".join ( "{} {}".format ( k + 1, p + 1 ) for k, p in enumerate ( sub["parent_index"] ) ) + "\n" } )
            skipped = [ ]
            mapped_cal = BL.HessianCalibration.mapped ( calibration, sub["parent_index"] ) if calibration is not None else None
            sem = SM.dyff_edits_from_hessian ( h, x_opt, [ ( a, b ) for a, b, o, ar in sub["bonds"] ], skipped, calibration = mapped_cal )
            return SS.map_term_edits_to_parent ( sem, sub["parent_index"], group ), len ( skipped ), len ( symbols )
        say ( "xTB Hessians of {} group(s) of atoms outside the library...".format ( len ( groups ) ) )
        sizes = [ ]
        with concurrent.futures.ThreadPoolExecutor ( max_workers = max ( 1, threads ) ) as pool:
            for k, ( mapped, n_skipped, n_sub ) in enumerate ( pool.map ( one, groups ), 1 ):
                edits["bond_terms"].update ( mapped["bond_terms"] ); edits["angle_terms"].update ( mapped["angle_terms"] )
                report["skipped"] += n_skipped; sizes.append ( n_sub )
                say ( "xTB Hessians {}/{}".format ( k, len ( groups ) ) )
        lap ( "hessian" )
        report["group_atoms"] = [ [ int ( a ) for a in g ] for g in groups ]
        report.update ( bonds = len ( edits["bond_terms"] ), angles = len ( edits["angle_terms"] ), groups = len ( groups ),
                        largest_subsystem = max ( sizes ) if sizes else 0,
                        calibrated = dict ( calibration.used ) if calibration is not None else { } )
        if torsions:
            say ( "torsion library..." )
            found, info = library_dihedrals ( z, bonds, formal, xyz, rest, charged = True )      # this model has charges
            edits["dihedral_bonds"].update ( found )
            report.update ( info )
            report["library_torsions"] = sorted ( found )
            if scan_torsions:
                # . 2026-10-21: the flexible bonds the library does not have -> xTB scan + fit
                targets = flexible_scan_targets ( z, bonds, formal, rest, set ( found ) )
                report["scan_targets"] = len ( targets )
                if targets:
                    fitted, failed, details = scan_flexible_torsions ( z, xyz, bonds, formal, groups, targets, charges, edits,
                                                              known_charges = known_charges, known_edits = known_edits,
                                                              shell = shell, solvent = solvent, multiplicity_of = multiplicity_of,
                                                              workers = max ( 1, threads ), say = say,
                                                              fragments = scan_fragments, fragment_gas = scan_gas,
                                                              log_dir = sub_log ( "torsions" ) )
                    edits["dihedral_bonds"].update ( fitted )
                    report["scanned"] = len ( fitted )
                    report["scan_failed"] = failed
                    report["scan_details"] = details
                    report["scan_settings"] = { "fragments": bool ( scan_fragments ), "fragment_shell": FRAGMENT_SHELL,
                                                "gas_phase": bool ( scan_fragments and scan_gas ), "step": 15.0,
                                                "min_improvement": SCAN_MIN_IMPROVEMENT, "poor_rmse": SCAN_POOR_RMSE }
                lap ( "torsions" )
    return { "charges": charges, "edits": edits, "report": report }


# . automatic torsion scans (scan_flexible_torsions): a fit is kept only when its RMSE against
#   the xTB profile is below this fraction of the DYFF rule's; above SCAN_POOR_RMSE (kcal/mol)
#   it is kept but reported (a profile that a torsion alone cannot describe, e.g. internal
#   hydrogen bonds that form and break along the relaxed scan)
SCAN_MIN_IMPROVEMENT = 0.9
SCAN_POOR_RMSE       = 2.0


def flexible_scan_targets ( z, bonds, formal, atoms, exclude = ( ) ):
    """ [ ( ( j, k ), ( i, j, k, l ) ) ] for the flexible bonds among `atoms` worth an xTB
        torsion scan: acyclic single bonds with heavy substituents on both sides, minus the
        rigid ones (in rings, multiple/aromatic, amide-like C(=O/S/N)-N with a planar N) and
        the symmetric terminal rotors (CH3, NH2, NH3+: a side whose other neighbours are two
        or more hydrogens). exclude: "j-k" keys already parametrized (torsion library).
        quad: the heaviest outer atoms (deterministic). """
    from collections import defaultdict
    members = set ( atoms )
    nb = defaultdict ( set ); order = { }; aromatic = set ( )
    for a, b, o, ar in bonds:
        nb[a].add ( b ); nb[b].add ( a ); order[( a, b )] = order[( b, a )] = int ( o )
        if ar: aromatic.add ( ( a, b ) ); aromatic.add ( ( b, a ) )
    def in_ring ( a, b ):
        seen, stack = { a }, [ a ]
        while stack:
            u = stack.pop ( )
            for v in nb[u]:
                if ( u, v ) in ( ( a, b ), ( b, a ) ): continue
                if v == b: return True
                if v not in seen: seen.add ( v ); stack.append ( v )
        return False
    def has_multiple ( a, to = None ):
        return any ( ( order[( a, v )] > 1 or ( a, v ) in aromatic ) and ( to is None or z[v] in to ) for v in nb[a] )
    def amide_like ( a, b ):
        for c, n in ( ( a, b ), ( b, a ) ):
            if z[n] == 7 and len ( nb[n] ) <= 3 and z[c] in ( 6, 15, 16 ) and has_multiple ( c, to = ( 7, 8, 16 ) ):
                return True
        return False
    def symmetric_rotor ( a, b ):
        outer = [ v for v in nb[a] if v != b ]
        return len ( outer ) >= 2 and all ( z[v] == 1 for v in outer )
    out = [ ]
    for a, b, o, ar in bonds:
        if a not in members or b not in members: continue
        j, k = min ( a, b ), max ( a, b )
        if "{}-{}".format ( j, k ) in exclude: continue
        if o != 1 or ar or in_ring ( j, k ) or amide_like ( j, k ): continue
        outer_j = [ v for v in nb[j] if v != k ]; outer_k = [ v for v in nb[k] if v != j ]
        if not outer_j or not outer_k: continue
        if symmetric_rotor ( j, k ) or symmetric_rotor ( k, j ): continue
        if all ( z[v] == 1 for v in outer_j + outer_k ): continue
        i = max ( outer_j, key = lambda v: ( z[v], -v ) ); l = max ( outer_k, key = lambda v: ( z[v], -v ) )
        out.append ( ( ( j, k ), ( i, j, k, l ) ) )
    return sorted ( out )


# . 2026-10-23: each rotatable bond is scanned on its own fragment: the two atoms of the
#   bond plus FRAGMENT_SHELL bonds around them, never cutting aromatic, multiple or small-
#   ring bonds nor metal bonds (subsystem.build_subsystem), capped with H. Smaller than the
#   capped subsystem of the whole ligand, so faster, and the far parts of the molecule can
#   no longer fold onto the scanned bond (internal hydrogen bonds that form and break along
#   a relaxed scan of the whole ligand).
FRAGMENT_SHELL = 2
# . 2026-10-23 (user decision): with fragments, the xTB scan runs in the GAS PHASE and the
#   MM side of the fit uses the fragment's own gas-phase CM5 charges (GFN1-xTB, caps
#   included, integer total). The fit's MM side is gas phase (NBModelFull, no screening),
#   so an ALPB scan put solvent-screened electrostatics on one side only: ATP mean fit
#   RMSE 4.5 (ALPB scans) -> 2.6 kcal/mol (gas scans + gas charges). The torsion is the
#   intrinsic one of the bond (as in GAFF / OpenFF fits); the medium enters through the
#   charges used in the simulation.


def torsion_fragment ( z, xyz, bonds, formal, bond, shell = FRAGMENT_SHELL ):
    """ The capped fragment for the torsion scan of `bond` (subsystem.build_subsystem dict). """
    from pdynamo import subsystem as SS
    return SS.build_subsystem ( z, xyz, bonds, list ( bond ), shell = shell, formal_charges = formal )


def scan_flexible_torsions ( z, xyz, bonds, formal, groups, targets, charges, edits, known_charges = None, known_edits = None,
                             shell = 2, solvent = "water", multiplicity_of = None, workers = 4, step = 15.0, say = None,
                             fragments = True, fragment_shell = FRAGMENT_SHELL, fragment_gas = True, log_dir = None ):
    """ Relaxed GFN2-xTB scans of the targets ( flexible_scan_targets ) on the capped
        subsystem of their group of rest atoms, and the self-consistent per-bond fit
        (pdynamo/dyff_bond_fit.py, as the Builder's DYFF Parameters "refine torsions")
        on top of the Seminario bonds/angles and the CM5 charges just computed.
        log_dir: per bond a folder torsion_<j+1>-<k+1> with the fragment and the original xtb outputs.
        Returns ( { "j-k": pairs } in system indices, [ "j-k: reason" ] for failures,
        [ { "bond", "before_rmse", "fit_rmse", "barrier_xtb", "accepted", ... } ] in kcal/mol;
        the rejected fits are in the details too, with accepted False ). """
    import contextlib
    import io
    from pScientific import PeriodicTable
    from pMolecule.NBModel import NBModelFull
    from pdynamo import subsystem as SS
    from pdynamo.dyff_edits import dyff_model, soften_polar_hydrogens
    from gui.windows.builder.dyff_parameters import _scan_and_fit
    say = say or ( lambda text: None )
    multiplicity_of = multiplicity_of or ( lambda atoms: 1 )
    all_charges = dict ( known_charges or { } ); all_charges.update ( charges )
    parent_edits = { kind: dict ( ( known_edits or { } ).get ( kind ) or { } ) for kind in ( "bond_terms", "angle_terms", "dihedral_bonds" ) }
    for kind in parent_edits: parent_edits[kind].update ( edits.get ( kind ) or { } )
    fitted, failed, details = { }, [ ], [ ]
    group_of = { a: g for g, group in enumerate ( groups ) for a in group }
    by_group = defaultdict ( list )
    for bond, quad in targets:
        if bond[0] in group_of: by_group[group_of[bond[0]]].append ( ( bond, quad ) )
    from pdynamo.dyff_bond_fit import STDOUT_LOCK as make_lock
    # . fragments=True: one scan system per bond (torsion_fragment), the fragments in parallel
    jobs = [ ]
    for g, items in sorted ( by_group.items ( ) ):
        group = groups[g]
        if fragments:
            for item in items:
                jobs.append ( ( group, [ item ], torsion_fragment ( z, xyz, bonds, formal, item[0], shell = fragment_shell ) ) )
        else:
            jobs.append ( ( group, items, SS.build_subsystem ( z, xyz, bonds, group, shell = shell, formal_charges = formal ) ) )
    def run_job ( job, inner_workers ):
        group, items, sub = job
        parent = sub["parent_index"]
        local = { p: k for k, p in enumerate ( parent ) if p >= 0 }
        local_targets = [ ]
        for bond, quad in items:
            if all ( a in local for a in quad ):
                q = tuple ( local[a] for a in quad )
                local_targets.append ( ( ( min ( q[1], q[2] ), max ( q[1], q[2] ) ), q ) )
            else:
                failed.append ( "{}-{}: dihedral cut by the subsystem".format ( *bond ) )
        if not local_targets: return
        symbols = [ PeriodicTable.Symbol ( x ) for x in sub["atomic_numbers"] ]
        pairs = [ ( a, b ) for a, b, o, ar in sub["bonds"] ]
        base = SS.map_term_edits_to_subsystem ( parent_edits, parent )
        def bond_dir ( b ):
            if not log_dir: return None
            j, k = sorted ( ( parent[b[0]], parent[b[1]] ) )
            return os.path.join ( log_dir, "torsion_{}-{}".format ( j + 1, k + 1 ) )
        job_dir = bond_dir ( local_targets[0][0] ) if fragments else ( os.path.join ( log_dir, "torsions_group_{}".format ( min ( group ) + 1 ) ) if log_dir else None )
        if job_dir:
            from pdynamo import param_logs
            param_logs.write_xyz ( os.path.join ( job_dir, "fragment.xyz" ), symbols, sub["coordinates"],
                                   "scan system: charge {} parent atoms (1-based, 0 = cap) {}".format (
                                       sub["charge"], " ".join ( str ( p + 1 ) for p in parent ) ) )
        sub_formal = [ formal[p] if p >= 0 else 0 for p in parent ]
        scan_solvent = solvent
        if fragments and fragment_gas:
            scan_solvent = None
            sub_charges = [ float ( q ) for q in xtb_cm5 ( symbols, sub["coordinates"], sub["charge"], solvent = None, threads = 1,
                                                           multiplicity = multiplicity_of ( set ( group ) ), log_dir = job_dir ) ]
        else:
            sub_charges = [ float ( all_charges.get ( p, 0.0 ) ) if p >= 0 else 0.0 for p in parent ]
        def make ( e, sub = sub, sub_formal = sub_formal, sub_charges = sub_charges ):
            # . the lock: redirect_stdout swaps the global sys.stdout, unsafe across threads
            with make_lock, contextlib.redirect_stdout ( io.StringIO ( ) ):
                system = SS.pdynamo_system ( sub["atomic_numbers"], sub["bonds"], sub["coordinates"], sub_formal )
                system.DefineMMModel ( dyff_model ( "dyff-1.0", e ), log = None )
                for k, q in enumerate ( sub_charges ): system.mmState.charges[k] = q
                soften_polar_hydrogens ( system )
                system.DefineNBModel ( NBModelFull.WithDefaults ( ) )
            return system
        def progress ( text, total = len ( targets ) ):
            say ( "torsion scans ({} bonds): {}".format ( total, text ) )
        results = _scan_and_fit ( symbols, sub["coordinates"], pairs, local_targets, make, base, sub["charge"],
                                  multiplicity_of ( set ( group ) ), scan_solvent, step, 1, inner_workers, progress,
                                  report_bond = lambda b, parent = parent: ( min ( parent[b[0]], parent[b[1]] ), max ( parent[b[0]], parent[b[1]] ) ),
                                  log_dir_of = ( bond_dir if fragments else ( lambda b: job_dir ) ) if log_dir else None )
        for r in results:
            key = "{}-{}".format ( *r["bond"] )
            if "pairs" in r:
                # . 2026-10-23: "after" = the re-relaxed DYFF profile with the fitted torsion, measured
                #   exactly like "before" (all points, unweighted); the fit's own weighted, best-offset
                #   RMSE made worse profiles look better (ATP P2-O8 6.11 -> 6.22 accepted as 4.06)
                before = r.get ( "before_rmse" )
                after = ( r.get ( "profiles" ) or { } ).get ( "after_rmse", r.get ( "fit_rmse" ) )
                accepted = not ( before is not None and after is not None and after > SCAN_MIN_IMPROVEMENT * before )
                entry = { "bond": key, "before_rmse": before, "fit_rmse": after, "barrier_xtb": r.get ( "barrier_xtb" ),
                          "poor": after is not None and after > SCAN_POOR_RMSE, "pairs": [ list ( p ) for p in r["pairs"] ],
                          "quad": [ parent[a] for a in r["quad"] ], "fragment_atoms": len ( sub["atomic_numbers"] ),
                          "fragment_charge": int ( sub["charge"] ), "scan_solvent": scan_solvent,
                          "dropped_periods": [ int ( x ) for x in ( r.get ( "dropped" ) or [ ] ) ],
                          "profiles": r.get ( "profiles" ), "accepted": accepted,
                          "log_dir": os.path.basename ( bond_dir ( r["quad"][1:3] ) ) if log_dir else None }
                details.append ( entry )
                if not accepted:
                    failed.append ( "{}: the fit does not improve on the DYFF rule ({:.2f} -> {:.2f} kcal/mol), rule kept".format ( key, before, after ) )
                    continue
                fitted[key] = entry["pairs"]
            else:
                failed.append ( "{}: {}".format ( key, r.get ( "error", "no fit" ) ) )
                details.append ( { "bond": key, "quad": [ parent[a] for a in r["quad"] ], "accepted": False, "error": r.get ( "error", "no fit" ),
                                   "fragment_atoms": len ( sub["atomic_numbers"] ),
                                   "log_dir": os.path.basename ( bond_dir ( r["quad"][1:3] ) ) if log_dir else None } )
    if fragments and len ( jobs ) > 1:
        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor ( max_workers = max ( 1, workers ) ) as pool:
            list ( pool.map ( lambda job: run_job ( job, 1 ), jobs ) )
    else:
        for job in jobs: run_job ( job, workers )
    return fitted, failed, details


def check_torsion ( z, xyz, bonds, formal, bond, edits, multiplicity = 1, step = 15.0, threads = 1, log_dir = None,
                    fragment_shell = FRAGMENT_SHELL ):
    """ 2026-10-23: how well a torsion ALREADY in `edits` (e.g. from the DFT torsion library)
        reproduces a relaxed GFN2-xTB scan, measured exactly as the automatic scans
        (scan_flexible_torsions): capped gas-phase fragment of the bond, its own gas-phase
        CM5 charges, DYFF relaxed at every point with the same restraint. Profiles of the
        DYFF rule (the bond's torsion edit removed) and of the given torsion.
        bond: ( j, k ) system indices; edits: dyff_parameter_edits-like dict (system indices).
        Returns { "bond", "quad", "fragment_atoms", "fragment_charge", "barrier_xtb", "rule_rmse",
        "torsion_rmse", "pairs", "profiles": { phi, xtb, dyff_rule, dyff_torsion } } in kcal/mol. """
    import contextlib
    import copy
    import io
    import numpy as np
    from pScientific import PeriodicTable
    from pMolecule.NBModel import NBModelFull
    from pdynamo import subsystem as SS
    from pdynamo import torsion_fit as TF
    from pdynamo import dyff_bond_fit as BF
    from pdynamo.dyff_edits import dyff_model, soften_polar_hydrogens
    j, k = min ( bond ), max ( bond )
    key = "{}-{}".format ( j, k )
    nb = defaultdict ( set )
    for a, b, *_ in bonds: nb[a].add ( b ); nb[b].add ( a )
    outer_j = [ v for v in nb[j] if v != k ]; outer_k = [ v for v in nb[k] if v != j ]
    if not outer_j or not outer_k: raise ValueError ( "bond {}-{} has no dihedral".format ( j + 1, k + 1 ) )
    quad = ( max ( outer_j, key = lambda v: ( z[v], -v ) ), j, k, max ( outer_k, key = lambda v: ( z[v], -v ) ) )
    sub = torsion_fragment ( z, xyz, bonds, formal, ( j, k ), shell = fragment_shell )
    parent = sub["parent_index"]
    local = { p: n for n, p in enumerate ( parent ) if p >= 0 }
    if not all ( a in local for a in quad ): raise ValueError ( "the dihedral is cut by the fragment" )
    q = tuple ( local[a] for a in quad )
    symbols = [ PeriodicTable.Symbol ( x ) for x in sub["atomic_numbers"] ]
    pairs_sub = [ ( a, b ) for a, b, o, ar in sub["bonds"] ]
    sub_formal = [ formal[p] if p >= 0 else 0 for p in parent ]
    if log_dir:
        from pdynamo import param_logs
        param_logs.write_xyz ( os.path.join ( log_dir, "fragment.xyz" ), symbols, sub["coordinates"],
                               "scan system: charge {} parent atoms (1-based, 0 = cap) {}".format (
                                   sub["charge"], " ".join ( str ( p + 1 ) for p in parent ) ) )
    sub_charges = [ float ( x ) for x in xtb_cm5 ( symbols, sub["coordinates"], sub["charge"], solvent = None, threads = threads,
                                                   multiplicity = multiplicity, log_dir = log_dir ) ]
    from pdynamo.dyff_bond_fit import STDOUT_LOCK
    def make ( e ):
        with STDOUT_LOCK, contextlib.redirect_stdout ( io.StringIO ( ) ):
            system = SS.pdynamo_system ( sub["atomic_numbers"], sub["bonds"], sub["coordinates"], sub_formal )
            system.DefineMMModel ( dyff_model ( "dyff-1.0", e ), log = None )
            for n, c in enumerate ( sub_charges ): system.mmState.charges[n] = c
            soften_polar_hydrogens ( system )
            system.DefineNBModel ( NBModelFull.WithDefaults ( ) )
        return system
    base = SS.map_term_edits_to_subsystem ( { kind: dict ( edits.get ( kind ) or { } ) for kind in ( "bond_terms", "angle_terms", "dihedral_bonds" ) }, parent )
    local_key = BF.bond_label ( q[1], q[2] )
    pairs = [ list ( x ) for x in ( edits.get ( "dihedral_bonds" ) or { } ).get ( key, [ ] ) ]
    with_torsion = copy.deepcopy ( base ); with_torsion.setdefault ( "dihedral_bonds", { } )[local_key] = pairs
    rule = copy.deepcopy ( base ); rule.setdefault ( "dihedral_bonds", { } ).pop ( local_key, None )
    runner = TF.XTBRunner ( gfn = 2, charge = sub["charge"], multiplicity = multiplicity, solvent = None, threads = threads, log_dir = log_dir )
    scan = TF.run_scan ( symbols, sub["coordinates"], pairs_sub, q, step = step, relaxed = True, runner = runner )
    e_qm = np.array ( scan["e_qm"], dtype = float )
    broken = TF.broken_bonds ( symbols, pairs_sub, scan["geometries"] )
    if broken:
        p, ( a, b ), d = broken[0]
        raise ValueError ( "the xTB scan changed the bonding (fragment atoms {}-{} at {:.2f} A)".format ( a + 1, b + 1, d ) )
    before = BF.profile_error ( make, rule, q, scan["geometries"], scan["phi"], e_qm )
    after = BF.profile_error ( make, with_torsion, q, scan["geometries"], scan["phi"], e_qm )
    return { "bond": key, "quad": [ int ( a ) for a in quad ], "fragment_atoms": len ( symbols ), "fragment_charge": int ( sub["charge"] ),
             "barrier_xtb": float ( ( e_qm.max ( ) - e_qm.min ( ) ) / 4.184 ), "pairs": pairs,
             "rule_rmse": before["rmse"] / 4.184, "torsion_rmse": after["rmse"] / 4.184,
             "log_dir": os.path.basename ( log_dir ) if log_dir else None,
             "profiles": { "phi": [ float ( x ) for x in scan["phi"] ], "xtb": [ float ( x ) for x in ( e_qm - e_qm.min ( ) ) / 4.184 ],
                           "dyff_rule": [ x / 4.184 for x in before["mm"] ], "dyff_torsion": [ x / 4.184 for x in after["mm"] ] } }


def library_dihedrals ( z, bonds, formal, xyz, atoms, charged = True ):
    """ Torsion corrections of the DFT library (pdynamo/bond_library.py) for the rotatable,
        acyclic single bonds whose atoms are all in `atoms`. charged: the table fitted with
        CM5 charges + polar-H LJ (the model with charges) or the one fitted without charges.
        Returns ( { "j-k": pairs }, report fields ). """
    from collections import defaultdict
    from pdynamo import bond_library as BL
    found_edits, report = { }, { "rotatable": 0, "torsions": 0, "levels": { } }
    library = ( BL.load_library ( charged = True ) or BL.load_library ( ) ) if charged else BL.load_library ( )
    if library is None: return found_edits, report
    members = set ( atoms )
    nb = defaultdict ( set )
    for a, b, *_ in bonds: nb[a].add ( b ); nb[b].add ( a )
    try:
        mol = BL.rdkit_molecule ( z, bonds, formal, xyz ) if BL.RDKIT_CODES else BL.code_molecule ( z, bonds, formal )
    except Exception as error:
        report["library_error"] = str ( error ).split ( "\n" )[0][:120]
        return found_edits, report
    targets = [ ]
    for a, b, order, aromatic in bonds:
        if a not in members or b not in members: continue
        if not BL.acyclic_single ( mol, a, b ): continue
        outer = [ u for u in nb[a] if u != b ] + [ u for u in nb[b] if u != a ]
        if not [ u for u in nb[a] if u != b ] or not [ u for u in nb[b] if u != a ]: continue
        if all ( z[u] == 1 for u in outer ): continue
        targets.append ( ( min ( a, b ), max ( a, b ) ) )
    report["rotatable"] = len ( targets )
    if targets:
        found = BL.library_torsions ( z, bonds, formal, targets, coordinates = xyz, library = library )
        for ( j, k ), ( level, pairs, count ) in found.items ( ):
            found_edits["{}-{}".format ( j, k )] = [ list ( p ) for p in pairs ]
            report["levels"][level] = report["levels"].get ( level, 0 ) + 1
        report["torsions"] = len ( found )
    return found_edits, report


def scale_charges ( charges, groups, scale ):
    """ CM5 charge scaling that keeps every group's total: q' = q_mean + scale * ( q - q_mean )
        within each group (a residue, a connected group of other atoms). Atoms in no
        group are left unchanged. Returns a new { atom: charge }. """
    out = dict ( charges )
    if scale == 1.0: return out
    for group in groups:
        group = [ a for a in group if a in charges ]
        if not group: continue
        mean = sum ( charges[a] for a in group ) / len ( group )
        for a in group: out[a] = mean + scale * ( charges[a] - mean )
    return out


def charge_groups ( matches, rest_atoms, neighbours, skip = ( ) ):
    """ The residues of the library matches and the connected groups of the other atoms
        (water and ions excluded through `skip`). """
    groups = [ sorted ( m["map"].values ( ) ) for m in matches ]
    members, seen = set ( rest_atoms ) - set ( skip ), set ( )
    for a in sorted ( members ):
        if a in seen: continue
        group, stack = [ ], [ a ]; seen.add ( a )
        while stack:
            u = stack.pop ( ); group.append ( u )
            for v in neighbours[u]:
                if v in members and v not in seen: seen.add ( v ); stack.append ( v )
        groups.append ( sorted ( group ) )
    return groups


def store_on_object ( vismol_object, system, edits, charges ):
    """ The library edits and charges on the object linked to the system (read by
        session.define_MMModel('DYFF')): per-term edits replaced, tied to the system's
        topology; every atom gets an mm_charge (0.0 where none). Returns the total charge. """
    from pdynamo.dyff_edits import term_signature
    current = getattr ( vismol_object, "dyff_parameter_edits", None ) or { }
    current = { k: v for k, v in current.items ( ) if k not in ( "bond_terms", "angle_terms", "dihedral_bonds", "term_signature" ) }
    current.update ( { "bond_terms": edits["bond_terms"], "angle_terms": edits["angle_terms"], "dihedral_bonds": edits["dihedral_bonds"],
                       "term_signature": term_signature ( len ( system.atoms ), list ( system.connectivity.bondIndices ) ) } )
    vismol_object.dyff_parameter_edits = current
    n = len ( system.atoms )
    values = [ float ( charges.get ( i, 0.0 ) ) for i in range ( n ) ]
    # . exact integer total (the library charges are stored with 5 decimals)
    residual = round ( sum ( values ) ) - sum ( values )
    values = [ v + residual / n for v in values ]
    for i, v in enumerate ( values ):
        atom = vismol_object.atoms.get ( i )
        if atom is not None: atom.mm_charge = v
    return sum ( values )

