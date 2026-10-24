#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: output folder of a DYFF-modified parametrization (2026-10-23)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Every run of "DYFF + residue library" / "Prepare DYFF System" writes a
#      folder next to the system (its working folder):
#
#        dyff_param_<system>_<date>/
#          parametrization.json   all the data: settings, per-atom charges (CM5,
#                                 RESP, final), Hessian bonds/angles, torsions
#                                 (library and fitted, with the scan profiles),
#                                 timings, warnings; read back by load() / apply()
#          report.html            the same data for people: tables, charge and
#                                 torsion plots (inline SVG), the original logs
#          logs/charges/cm5/      xtb CM5 run (subsystem of the atoms outside the library)
#          logs/charges/resp_N/   xtb --esp or ORCA (+ orca_vpot) run of the molecule
#                                 whose first atom is N
#          logs/hessian/group_N/  xtb --ohess of the subsystem of group N
#          logs/torsions/torsion_J-K/  fragment, CM5, relaxed scan (xtbscan.log is a
#                                 multi-frame xyz) and the single points
#
#      Atom numbers in file and folder names and in the report are 1-based (as
#      in the EasyHybrid interface); the keys inside parametrization.json are
#      0-based system indices (as in dyff_parameter_edits).
#
#      JSON, not pickle: a pickle ties the file to the classes and numpy version
#      that wrote it (see the .easy numpy compatibility shim); JSON stays readable
#      by any version and by other programs.
#

import datetime
import html
import json
import math
import os

FORMAT  = "easyhybrid-dyff-parametrization"
VERSION = 1
KCAL    = 4.184
CMAX_EMBED = 400000           # bytes: larger logs (the 24 single points of a scan) are linked, not embedded


#=============================================================================
#  Folder
#=============================================================================
def new_folder ( base, name ):
    """ <base>/dyff_param_<name>_<date-time> (created). """
    safe = "".join ( c if c.isalnum ( ) or c in "-_." else "_" for c in ( name or "system" ) )[:60]
    stamp = datetime.datetime.now ( ).strftime ( "%Y%m%d-%H%M%S" )
    folder = os.path.join ( base or os.getcwd ( ), "dyff_param_{}_{}".format ( safe, stamp ) )
    os.makedirs ( os.path.join ( folder, "logs" ), exist_ok = True )
    return folder


def _version ( ):
    try:
        from _version import EASYHYBRID_VERSION
        return str ( EASYHYBRID_VERSION )
    except Exception:
        return "unknown"


def _atom_label ( atom, index ):
    path = getattr ( atom, "path", None ) or ""
    parts = [ p for p in path.split ( ":" ) if p ]
    name = parts[-1] if parts else getattr ( atom, "label", "" ) or "?"
    residue = parts[-2] if len ( parts ) >= 2 else ""
    return name, residue


#=============================================================================
#  Writing
#=============================================================================
def collect ( system, plan, rest, settings, charges, edits, messages = None, total_seconds = None ):
    """ The parametrization.json dictionary. plan: residue_library.plan; rest: the
        rest_parameters result (or None); charges / edits: the FINAL ones applied
        (library + rest, after scaling); settings: the dialog options. """
    from pScientific import PeriodicTable
    from pdynamo.dyff_edits import term_signature
    n = len ( system.atoms )
    report = dict ( ( rest or { } ).get ( "report" ) or { } )
    rest_atoms = sorted ( plan.get ( "rest_atoms" ) or [ ] )
    cm5 = { int ( k ): v for k, v in ( report.pop ( "cm5_charges", None ) or { } ).items ( ) }
    resp, esp = { }, { }
    for g in report.get ( "ligand_charges" ) or [ ]:
        if g.get ( "status", "" ).startswith ( "RESP" ):
            for a, q in zip ( g["members"], g.get ( "esp_charges" ) or [ ] ): esp[a] = q
            for a in g["members"]: resp[a] = ( rest or { } ).get ( "charges", { } ).get ( a )
    atoms = [ ]
    for a in rest_atoms:
        atom = system.atoms[a]
        name, residue = _atom_label ( atom, a )
        atoms.append ( { "index": a, "name": name, "residue": residue, "element": PeriodicTable.Symbol ( atom.atomicNumber ),
                         "formal_charge": int ( getattr ( atom, "formalCharge", 0 ) or 0 ),
                         "xyz": [ round ( float ( system.coordinates3[a, k] ), 5 ) for k in range ( 3 ) ],
                         "cm5": cm5.get ( a ), "resp": resp.get ( a ), "esp": esp.get ( a ), "charge": charges.get ( a ) } )
    rest_set = set ( rest_atoms )
    rest_edits = { kind: { key: value for key, value in ( edits.get ( kind ) or { } ).items ( )
                           if all ( int ( x ) in rest_set for x in key.split ( "-" ) ) }
                   for kind in ( "bond_terms", "angle_terms", "dihedral_bonds" ) }
    bonds = [ ( b.node1.index, b.node2.index ) for b in system.connectivity.bonds ]
    s = plan.get ( "summary" ) or { }
    return { "format": FORMAT, "format_version": VERSION, "easyhybrid_version": _version ( ),
             "created": datetime.datetime.now ( ).isoformat ( timespec = "seconds" ),
             "system": { "name": getattr ( system, "label", "" ) or "", "atoms": n,
                         "elements": "".join ( PeriodicTable.Symbol ( a.atomicNumber ) + " " for a in system.atoms ).split ( ),
                         "term_signature": term_signature ( n, bonds ),
                         "total_charge": round ( sum ( float ( charges.get ( i, 0.0 ) ) for i in range ( n ) ), 6 ) + 0.0 },
             "settings": settings, "library": { "residues": s.get ( "residues" ), "matched": s.get ( "matched" ),
                                                 "templates": s.get ( "templates" ), "unmatched": s.get ( "unmatched" ),
                                                 "waters": s.get ( "waters" ), "ions": s.get ( "ions" ) },
             "rest_atoms": rest_atoms, "atoms": atoms,
             "charges": { str ( a ): charges.get ( a ) for a in rest_atoms },
             "rest_charges": { str ( a ): q for a, q in ( ( rest or { } ).get ( "charges" ) or { } ).items ( ) },
             "edits": rest_edits, "report": report, "messages": messages or [ ],
             "total_seconds": total_seconds }


def write ( folder, data ):
    """ parametrization.json + report.html in `folder`. Returns the report path. """
    with open ( os.path.join ( folder, "parametrization.json" ), "w" ) as handle:
        json.dump ( data, handle, indent = 1, default = _json_default )
    path = os.path.join ( folder, "report.html" )
    with open ( path, "w" ) as handle:
        handle.write ( render_html ( data, folder ) )
    return path


def _json_default ( value ):
    try:
        import numpy as np
        if isinstance ( value, np.generic ): return value.item ( )
        if isinstance ( value, np.ndarray ): return value.tolist ( )
    except Exception:
        pass
    if isinstance ( value, ( set, tuple ) ): return list ( value )
    return str ( value )


#=============================================================================
#  Reading back
#=============================================================================
def load ( folder_or_file ):
    path = folder_or_file
    if os.path.isdir ( path ): path = os.path.join ( path, "parametrization.json" )
    with open ( path ) as handle: data = json.load ( handle )
    if data.get ( "format" ) != FORMAT:
        raise ValueError ( "{} is not an EasyHybrid DYFF parametrization file.".format ( path ) )
    data["_folder"] = os.path.dirname ( os.path.abspath ( path ) )
    return data


def check_system ( data, system ):
    """ None when the stored parameters fit this system (same atoms, elements and
        bonds), else the reason. """
    from pScientific import PeriodicTable
    from pdynamo.dyff_edits import term_signature
    n = len ( system.atoms )
    if n != data["system"]["atoms"]:
        return "the system has {} atoms, the parametrization {}".format ( n, data["system"]["atoms"] )
    elements = [ PeriodicTable.Symbol ( a.atomicNumber ) for a in system.atoms ]
    if elements != data["system"]["elements"]:
        first = next ( i for i, ( a, b ) in enumerate ( zip ( elements, data["system"]["elements"] ) ) if a != b )
        return "different element at atom {} ({} here, {} in the parametrization)".format ( first + 1, elements[first], data["system"]["elements"][first] )
    bonds = [ ( b.node1.index, b.node2.index ) for b in system.connectivity.bonds ]
    if term_signature ( n, bonds ) != data["system"]["term_signature"]:
        return "the bonds of the system differ from the parametrized ones"
    return None


def rest_result ( data ):
    """ The stored parameters in the form of residue_library.rest_parameters'
        result, so the caller can apply them exactly as a fresh calculation. """
    # . rest_charges: before the optional CM5 scaling, which the caller applies again
    return { "charges": { int ( k ): v for k, v in ( data.get ( "rest_charges" ) or data["charges"] ).items ( ) if v is not None },
             "edits": { kind: dict ( data["edits"].get ( kind ) or { } ) for kind in ( "bond_terms", "angle_terms", "dihedral_bonds" ) },
             "report": data.get ( "report" ) or { } }


#=============================================================================
#  SVG plots (no matplotlib in the GUI's Python)
#=============================================================================
def _ticks ( lo, hi, n = 5 ):
    if hi <= lo: hi = lo + 1.0
    raw = ( hi - lo ) / n
    mag = 10 ** math.floor ( math.log10 ( raw ) )
    step = min ( ( m * mag for m in ( 1, 2, 2.5, 5, 10 ) ), key = lambda s: abs ( s - raw ) )
    start = math.floor ( lo / step ) * step
    out, v = [ ], start
    while v <= hi + 1e-9:
        out.append ( round ( v, 10 ) ); v += step
    return out


def svg_lines ( series, xlabel, ylabel, xticks = None, width = 520, height = 300, ylim = None ):
    """ series: [ ( label, xs, ys, css class, dashed, markers ) ]. """
    ml, mr, mt, mb = 52, 12, 12, 40
    xs_all = [ x for s in series for x in s[1] ]; ys_all = [ y for s in series for y in s[2] if y is not None ]
    if not xs_all or not ys_all: return ""
    x0, x1 = ( min ( xticks ), max ( xticks ) ) if xticks else ( min ( xs_all ), max ( xs_all ) )
    y0, y1 = ylim or ( min ( 0.0, min ( ys_all ) ), max ( ys_all ) * 1.08 + 1e-6 )
    yt = _ticks ( y0, y1 ); y0, y1 = min ( y0, yt[0] ), max ( y1, yt[-1] )
    X = lambda x: ml + ( x - x0 ) / ( x1 - x0 ) * ( width - ml - mr )
    Y = lambda y: mt + ( y1 - y ) / ( y1 - y0 ) * ( height - mt - mb )
    out = [ '<svg class="plot" viewBox="0 0 {} {}" role="img">'.format ( width, height ) ]
    for t in yt:
        out.append ( '<line class="grid" x1="{:.1f}" x2="{:.1f}" y1="{:.1f}" y2="{:.1f}"/>'.format ( ml, width - mr, Y ( t ), Y ( t ) ) )
        out.append ( '<text class="tick" x="{:.1f}" y="{:.1f}" text-anchor="end">{:g}</text>'.format ( ml - 5, Y ( t ) + 4, t ) )
    for t in ( xticks or _ticks ( x0, x1 ) ):
        out.append ( '<text class="tick" x="{:.1f}" y="{:.1f}" text-anchor="middle">{:g}</text>'.format ( X ( t ), height - mb + 15, t ) )
    out.append ( '<line class="axis" x1="{0}" x2="{0}" y1="{1}" y2="{2}"/>'.format ( ml, mt, height - mb ) )
    out.append ( '<line class="axis" x1="{0}" x2="{1}" y1="{2}" y2="{2}"/>'.format ( ml, width - mr, height - mb ) )
    out.append ( '<text class="label" x="{:.1f}" y="{}" text-anchor="middle">{}</text>'.format ( ( ml + width - mr ) / 2, height - 6, html.escape ( xlabel ) ) )
    out.append ( '<text class="label" transform="translate(13,{:.1f}) rotate(-90)" text-anchor="middle">{}</text>'.format ( ( mt + height - mb ) / 2, html.escape ( ylabel ) ) )
    for label, xs, ys, css, dashed, markers in series:
        pts = [ ( X ( x ), Y ( y ) ) for x, y in zip ( xs, ys ) if y is not None ]
        out.append ( '<polyline class="s {}{}" points="{}"/>'.format ( css, " dash" if dashed else "", " ".join ( "{:.1f},{:.1f}".format ( *p ) for p in pts ) ) )
        if markers:
            out += [ '<circle class="m {}" cx="{:.1f}" cy="{:.1f}" r="2.6"/>'.format ( css, *p ) for p in pts ]
    out.append ( "</svg>" )
    return "".join ( out )


def svg_bars ( labels, groups, ylabel, width = None, height = 260 ):
    """ groups: [ ( label, values, css class ) ] side by side for every label. """
    ml, mr, mt, mb = 52, 10, 12, 58
    n, g = len ( labels ), len ( groups )
    if not n: return ""
    width = width or max ( 360, ml + mr + n * ( 10 + 9 * g ) )
    vals = [ v for _, vs, _ in groups for v in vs if v is not None ] + [ 0.0 ]
    yt = _ticks ( min ( vals ), max ( vals ) ); y0, y1 = yt[0], yt[-1]
    Y = lambda y: mt + ( y1 - y ) / ( y1 - y0 ) * ( height - mt - mb )
    slot = ( width - ml - mr ) / n; bw = slot * 0.8 / g
    out = [ '<svg class="plot" viewBox="0 0 {:.0f} {}" role="img">'.format ( width, height ) ]
    for t in yt:
        out.append ( '<line class="grid" x1="{}" x2="{:.1f}" y1="{:.1f}" y2="{:.1f}"/>'.format ( ml, width - mr, Y ( t ), Y ( t ) ) )
        out.append ( '<text class="tick" x="{}" y="{:.1f}" text-anchor="end">{:g}</text>'.format ( ml - 5, Y ( t ) + 4, t ) )
    out.append ( '<line class="axis" x1="{0}" x2="{1:.1f}" y1="{2:.1f}" y2="{2:.1f}"/>'.format ( ml, width - mr, Y ( 0.0 ) ) )
    for i, lab in enumerate ( labels ):
        for k, ( _, vs, css ) in enumerate ( groups ):
            v = vs[i]
            if v is None: continue
            x = ml + i * slot + slot * 0.1 + k * bw
            out.append ( '<rect class="b {}" x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}"><title>{}: {:.3f}</title></rect>'.format (
                css, x, min ( Y ( v ), Y ( 0 ) ), max ( bw - 0.6, 0.6 ), abs ( Y ( v ) - Y ( 0 ) ), html.escape ( lab ), v ) )
        out.append ( '<text class="tick" transform="translate({:.1f},{}) rotate(-60)" text-anchor="end">{}</text>'.format (
            ml + ( i + 0.5 ) * slot, height - mb + 10, html.escape ( lab ) ) )
    out.append ( '<text class="label" transform="translate(13,{:.1f}) rotate(-90)" text-anchor="middle">{}</text>'.format ( ( mt + height - mb ) / 2, html.escape ( ylabel ) ) )
    out.append ( "</svg>" )
    return "".join ( out )


#=============================================================================
#  HTML
#=============================================================================
CSS = """
:root{--bg:#ffffff;--fg:#1d1d1f;--muted:#5a6270;--line:#dde2ea;--head:#f1f4f8;--accent:#1565c0;--box:#f4f7fb;
--warn:#ef6c00;--warnbg:#fff5ea;--good:#2e7d32;--bad:#c62828;--rule:#757575;--code:#eef1f5;--s1:#1d1d1f;--s2:#757575;--s3:#1565c0;--s4:#c62828;--s5:#ef6c00}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#14171c;--fg:#e6e9ee;--muted:#9aa3b1;--line:#2c323b;--head:#1d222a;
--accent:#6aa9f0;--box:#1a2029;--warn:#ffa25a;--warnbg:#2a2018;--good:#7cc77f;--bad:#ff7b72;--rule:#9aa3b1;--code:#232931;--s1:#e6e9ee;--s2:#9aa3b1;--s3:#6aa9f0;--s4:#ff7b72;--s5:#ffa25a;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#14171c;--fg:#e6e9ee;--muted:#9aa3b1;--line:#2c323b;--head:#1d222a;--accent:#6aa9f0;--box:#1a2029;--warn:#ffa25a;--warnbg:#2a2018;
--good:#7cc77f;--bad:#ff7b72;--rule:#9aa3b1;--code:#232931;--s1:#e6e9ee;--s2:#9aa3b1;--s3:#6aa9f0;--s4:#ff7b72;--s5:#ffa25a;color-scheme:dark}
body{background:var(--bg);color:var(--fg);font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;font-size:15px;line-height:1.5;margin:0;
padding-block:24px 60px;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
main{max-width:1060px;margin:0 auto}
h1{font-size:1.6em;margin:0 0 4px}h2{font-size:1.22em;border-bottom:2px solid var(--accent);padding-bottom:4px;margin-top:2em}
h3{font-size:1.05em;margin:1.4em 0 .4em}.meta,.small{color:var(--muted);font-size:.88em}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px;margin:14px 0}
.card{background:var(--box);border-left:4px solid var(--accent);padding:8px 12px}.card b{display:block;font-size:1.3em}
.box{background:var(--box);border-left:4px solid var(--accent);padding:8px 14px;margin:12px 0}
.box.warn{border-left-color:var(--warn);background:var(--warnbg)}
.tw{overflow-x:auto;margin:10px 0}table{border-collapse:collapse;font-size:.88em;min-width:50%}
th,td{border:1px solid var(--line);padding:3px 7px;text-align:right;white-space:nowrap}th{background:var(--head)}
td.l,th.l{text-align:left}.ok{color:var(--good);font-weight:600}.no{color:var(--bad);font-weight:600}
code,pre{font-family:ui-monospace,"SFMono-Regular",Menlo,Consolas,monospace;background:var(--code);border-radius:3px}
code{padding:1px 4px;font-size:.88em}pre{padding:8px;font-size:.78em;overflow:auto;max-height:420px}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,480px),1fr));gap:18px}
.torsion{border:1px solid var(--line);border-radius:6px;padding:10px 14px;min-width:0}
svg.plot{width:100%;height:auto;max-width:100%}
svg .grid{stroke:var(--line);stroke-width:1}svg .axis{stroke:var(--muted);stroke-width:1}
svg .tick{fill:var(--muted);font-size:10px}svg .label{fill:var(--fg);font-size:11px}
svg .s{fill:none;stroke-width:2}svg .dash{stroke-dasharray:5 4}
.c1{stroke:var(--s1);fill:var(--s1)}.c2{stroke:var(--s2);fill:var(--s2)}.c3{stroke:var(--s3);fill:var(--s3)}.c4{stroke:var(--s4);fill:var(--s4)}.c5{stroke:var(--s5);fill:var(--s5)}
svg .s.c1,svg .s.c2,svg .s.c3,svg .s.c4,svg .s.c5{fill:none}
.legend{display:flex;flex-wrap:wrap;gap:12px;font-size:.85em;color:var(--muted)}.legend i{display:inline-block;width:18px;height:3px;vertical-align:middle;margin-right:4px}
.k1{background:var(--s1)}.k2{background:var(--s2)}.k3{background:var(--s3)}.k4{background:var(--s4)}.k5{background:var(--s5)}
details{margin:4px 0}summary{cursor:pointer}
"""


def _f ( v, fmt = "{:.3f}" ):
    return "" if v is None else fmt.format ( v )


def _alabel ( data, index ):
    atom = data["_by_index"].get ( index )
    if atom is None: return str ( index + 1 )
    return "{} ({})".format ( atom["name"], index + 1 )


def render_html ( data, folder ):
    data = dict ( data ); data["_by_index"] = { a["index"]: a for a in data["atoms"] }
    rep = data.get ( "report" ) or { }
    st = data.get ( "settings" ) or { }
    name = data["system"]["name"] or "system"
    E = html.escape
    out = [ "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
            "<title>DYFF parametrization: {}</title><style>{}</style></head><body><main>".format ( E ( name ), CSS ) ]
    out.append ( "<h1>DYFF parametrization: {}</h1>".format ( E ( name ) ) )
    out.append ( '<p class="meta">EasyHybrid {} · {} · folder <code>{}</code></p>'.format ( E ( data["easyhybrid_version"] ), E ( data["created"] ), E ( os.path.basename ( folder ) ) ) )
    # . summary
    details = rep.get ( "scan_details" ) or [ ]
    accepted = [ d for d in details if d.get ( "accepted" ) ]
    cards = [ ( "atoms in the system", data["system"]["atoms"] ), ( "atoms outside the library", len ( data["rest_atoms"] ) ),
              ( "total charge", "{:+.3f}".format ( round ( data["system"]["total_charge"], 3 ) + 0.0 ) ),
              ( "Hessian bonds / angles", "{} / {}".format ( rep.get ( "bonds", 0 ), rep.get ( "angles", 0 ) ) ),
              ( "library torsions", rep.get ( "torsions", 0 ) ),
              ( "scanned torsions fitted", "{} of {}".format ( len ( accepted ), rep.get ( "scan_targets", 0 ) ) ) ]
    if data.get ( "total_seconds" ) is not None: cards.append ( ( "wall time", "{:.0f} s".format ( data["total_seconds"] ) ) )
    out.append ( '<div class="cards">' + "".join ( '<div class="card">{}<b>{}</b></div>'.format ( E ( k ), E ( str ( v ) ) ) for k, v in cards ) + "</div>" )
    if data.get ( "messages" ):
        out.append ( '<div class="box warn">' + "".join ( "<p>{}</p>".format ( E ( m ) ) for m in data["messages"] ) + "</div>" )
    # . settings
    out.append ( '<h2 id="settings">Settings</h2><div class="tw"><table>' )
    lib = data.get ( "library" ) or { }
    rows = [ ( "residue library", "{} of {} residues matched{}".format ( lib.get ( "matched" ), lib.get ( "residues" ),
                                    "; unmatched: " + ", ".join ( lib.get ( "unmatched" ) or [ ] ) if lib.get ( "unmatched" ) else "" ) ) ]
    rows += [ ( k.replace ( "_", " " ), json.dumps ( v ) if isinstance ( v, ( dict, list ) ) else str ( v ) ) for k, v in st.items ( ) ]
    for k, v in ( rep.get ( "charge_settings" ) or { } ).items ( ): rows.append ( ( "charges: " + k.replace ( "_", " " ), str ( v ) ) )
    for k, v in ( rep.get ( "scan_settings" ) or { } ).items ( ): rows.append ( ( "torsion scans: " + k.replace ( "_", " " ), str ( v ) ) )
    if rep.get ( "timings" ):
        rows.append ( ( "wall time per step (s)", ", ".join ( "{} {}".format ( k, v ) for k, v in rep["timings"].items ( ) ) ) )
    out += [ '<tr><th class="l">{}</th><td class="l">{}</td></tr>'.format ( E ( a ), E ( b ) ) for a, b in rows ]
    out.append ( "</table></div>" )
    out.append ( _charges_html ( data, rep ) )
    out.append ( _bonded_html ( data, rep ) )
    out.append ( _torsions_html ( data, rep ) )
    out.append ( _logs_html ( folder ) )
    out.append ( '<h2 id="data">Data files</h2><ul><li><code>parametrization.json</code>: everything in this report '
                 '(format {} v{}). <b>System &gt; Force Field &gt; Load DYFF Parametrization...</b> applies it again to the same system '
                 'without recomputing; terminal: <code>dyffload folder=...</code>. Keys are 0-based system indices.</li>'
                 '<li><code>logs/</code>: the original inputs and outputs of xtb and ORCA. <code>xtbscan.log</code> files are multi-frame '
                 'xyz trajectories of the relaxed scans.</li></ul>'.format ( FORMAT, VERSION ) )
    out.append ( "</main></body></html>" )
    return "\n".join ( out )


def _charges_html ( data, rep ):
    E = html.escape
    out = [ '<h2 id="charges">Charges</h2>' ]
    groups = rep.get ( "ligand_charges" ) or [ ]
    model = ( rep.get ( "charge_settings" ) or { } ).get ( "ligand_charges", "cm5" )
    out.append ( "<p>Atoms outside the residue library: {} on a capped subsystem, integer total per connected group.".format (
        E ( ( rep.get ( "charge_settings" ) or { } ).get ( "cm5", "GFN1-xTB CM5" ) ) ) +
        ( " Ligand charge model: <b>{}</b>.".format ( E ( model ) ) if model != "cm5" else "" ) + "</p>" )
    if groups:
        out.append ( '<div class="tw"><table><tr><th class="l">molecule (first atom)</th><th>atoms</th><th class="l">charges</th><th>engine</th>'
                     '<th>charge</th><th>mult.</th><th>RRMS</th><th>dipole QM / RESP (D)</th><th>points</th><th>time (s)</th><th class="l">logs</th></tr>' )
        for g in groups:
            first = g["members"][0] if g.get ( "members" ) else None
            out.append ( '<tr><td class="l">{}</td><td>{}</td><td class="l">{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td class="l">{}</td></tr>'.format (
                E ( _alabel ( data, first ) ) if first is not None else "", g["atoms"], E ( g["status"] ), E ( str ( g.get ( "engine" ) or "" ) ),
                g.get ( "charge", "" ), g.get ( "multiplicity", "" ), _f ( g.get ( "rrms" ), "{:.4f}" ),
                ( "{} / {}".format ( _f ( g.get ( "dipole_qm" ), "{:.2f}" ), _f ( g.get ( "dipole_resp" ), "{:.2f}" ) ) if g.get ( "dipole_resp" ) is not None else "" ),
                g.get ( "points" ) or "", g.get ( "seconds", "" ),
                '<a href="logs/charges/{0}/">{0}</a>'.format ( E ( g["log_dir"] ) ) if g.get ( "log_dir" ) else "" ) )
        out.append ( "</table></div>" )
    atoms = data["atoms"]
    if atoms and len ( atoms ) <= 400:
        has_resp = any ( a.get ( "resp" ) is not None for a in atoms )
        labels = [ "{}{}".format ( a["name"], a["index"] + 1 ) for a in atoms ]
        series = [ ( "CM5", [ a.get ( "cm5" ) for a in atoms ], "c2" ) ]
        if has_resp: series.append ( ( "RESP-2", [ a.get ( "resp" ) for a in atoms ], "c3" ) )
        else: series = [ ( "final", [ a.get ( "charge" ) for a in atoms ], "c3" ) ]
        out.append ( '<div class="legend">' + "".join ( '<span><i class="k{}"></i>{}</span>'.format ( css[1], E ( l ) ) for l, _, css in series ) + "</div>" )
        out.append ( '<div class="tw">' + svg_bars ( labels, series, "charge (e)", width = max ( 480, 70 + len ( atoms ) * ( 8 + 8 * len ( series ) ) ) ) + "</div>" )
    out.append ( '<details><summary>Per-atom table ({} atoms)</summary><div class="tw"><table><tr><th>#</th><th class="l">atom</th><th class="l">residue</th>'
                 '<th class="l">el.</th><th>formal</th><th>CM5</th><th>ESP</th><th>RESP-2</th><th>final</th></tr>'.format ( len ( atoms ) ) )
    for a in atoms:
        out.append ( '<tr><td>{}</td><td class="l">{}</td><td class="l">{}</td><td class="l">{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td><b>{}</b></td></tr>'.format (
            a["index"] + 1, E ( a["name"] ), E ( a["residue"] ), a["element"], a["formal_charge"] or "", _f ( a.get ( "cm5" ), "{:.4f}" ),
            _f ( a.get ( "esp" ), "{:.4f}" ), _f ( a.get ( "resp" ), "{:.4f}" ), _f ( a.get ( "charge" ), "{:.4f}" ) ) )
    out.append ( "</table></div></details>" )
    return "\n".join ( out )


def _bonded_html ( data, rep ):
    E = html.escape
    edits = data.get ( "edits" ) or { }
    bonds, angles = edits.get ( "bond_terms" ) or { }, edits.get ( "angle_terms" ) or { }
    out = [ '<h2 id="bonded">Bonds and angles (modified Seminario, xTB Hessian)</h2>' ]
    out.append ( "<p>{} bond and {} angle terms from the Hessian of {} capped group(s) (largest subsystem {} atoms); {} term(s) left to the DYFF rule (non-physical constant).{}</p>".format (
        len ( bonds ), len ( angles ), rep.get ( "groups", 0 ), rep.get ( "largest_subsystem", 0 ), rep.get ( "skipped", 0 ),
        " DFT calibration by environment: " + E ( json.dumps ( rep["calibrated"] ) ) if rep.get ( "calibrated" ) else "" ) )
    if bonds:
        out.append ( '<details><summary>Bond terms ({})</summary><div class="tw"><table><tr><th class="l">bond</th><th>b0 (Å)</th><th>k (kJ/mol/Å²), E = k(b−b0)²</th></tr>'.format ( len ( bonds ) ) )
        for key, ( b0, k ) in sorted ( bonds.items ( ), key = lambda t: [ int ( x ) for x in t[0].split ( "-" ) ] ):
            out.append ( '<tr><td class="l">{}</td><td>{:.4f}</td><td>{:.1f}</td></tr>'.format ( E ( " – ".join ( _alabel ( data, int ( x ) ) for x in key.split ( "-" ) ) ), b0, k ) )
        out.append ( "</table></div></details>" )
    if angles:
        out.append ( '<details><summary>Angle terms ({})</summary><div class="tw"><table><tr><th class="l">angle</th><th>θ0 (°)</th><th>k (kJ/mol/rad²)</th></tr>'.format ( len ( angles ) ) )
        for key, pairs in sorted ( angles.items ( ), key = lambda t: [ int ( x ) for x in t[0].split ( "-" ) ] ):
            theta0, k = _angle_from_pairs ( pairs )
            out.append ( '<tr><td class="l">{}</td><td>{}</td><td>{}</td></tr>'.format ( E ( " – ".join ( _alabel ( data, int ( x ) ) for x in key.split ( "-" ) ) ),
                         _f ( theta0, "{:.2f}" ), _f ( k, "{:.1f}" ) ) )
        out.append ( "</table></div></details>" )
    return "\n".join ( out )


def _angle_from_pairs ( pairs ):
    """ ( theta0 deg, k ) back from seminario.dyff_angle_pairs. """
    c = { int ( n ): float ( v ) for v, n in pairs }
    if 2 not in c: return 180.0, c.get ( 1 )
    cos0 = max ( -1.0, min ( 1.0, -c.get ( 1, 0.0 ) / ( 4.0 * c[2] ) ) ) if c[2] else None
    if cos0 is None: return None, None
    theta0 = math.acos ( cos0 )
    return math.degrees ( theta0 ), 4.0 * c[2] * math.sin ( theta0 ) ** 2


def _torsions_html ( data, rep ):
    E = html.escape
    out = [ '<h2 id="torsions">Torsions</h2>' ]
    dih = ( data.get ( "edits" ) or { } ).get ( "dihedral_bonds" ) or { }
    library = rep.get ( "library_torsions" ) or [ ]
    out.append ( "<p>{} rotatable bond(s) among the atoms outside the library; {} from the DFT torsion library ({}).</p>".format (
        rep.get ( "rotatable", 0 ), rep.get ( "torsions", 0 ), E ( ", ".join ( "{} {}".format ( v, k ) for k, v in ( rep.get ( "levels" ) or { } ).items ( ) ) ) or "none" ) )
    if library:
        out.append ( '<details><summary>Library torsions ({})</summary><div class="tw"><table><tr><th class="l">bond</th><th class="l">series E = Σ c cos(nφ) (c kJ/mol, n)</th></tr>'.format ( len ( library ) ) )
        for key in library:
            out.append ( '<tr><td class="l">{}</td><td class="l">{}</td></tr>'.format ( E ( " – ".join ( _alabel ( data, int ( x ) ) for x in key.split ( "-" ) ) ),
                         E ( "  ".join ( "{:+.3f}×{}".format ( c, int ( n ) ) for c, n in dih.get ( key, [ ] ) ) ) ) )
        out.append ( "</table></div></details>" )
    checks = { k: c for k, c in ( rep.get ( "torsion_checks" ) or { } ).items ( ) if c.get ( "profiles" ) }
    if checks:
        # . 2026-10-23: library torsions checked against an xTB scan (Parametrization window)
        out.append ( "<h3>Library torsions checked against xTB scans</h3><p>Same relaxed GFN2-xTB fragment scan as the automatic fits. "
                     "The library torsions were fitted to DFT, so this is a check of consistency with xTB, not a fit.</p>" )
        out.append ( '<div class="legend"><span><i class="k1"></i>GFN2-xTB</span><span><i class="k2"></i>DYFF rule</span>'
                     '<span><i class="k5"></i>DYFF + library torsion</span></div><div class="grid2">' )
        for key, c in sorted ( checks.items ( ), key = lambda t: [ int ( x ) for x in t[0].split ( "-" ) ] ):
            j, k = ( int ( x ) for x in key.split ( "-" ) )
            p = c["profiles"]
            phi = [ ( ( x + 180.0 ) % 360.0 ) - 180.0 for x in p["phi"] ]
            order = sorted ( range ( len ( phi ) ), key = lambda i: phi[i] )
            xs = [ phi[i] for i in order ]
            pick = lambda name: [ p[name][i] for i in order ]
            out.append ( '<div class="torsion"><h3>{} – {} (library)</h3>'.format ( E ( _alabel ( data, j ) ), E ( _alabel ( data, k ) ) ) )
            out.append ( svg_lines ( [ ( "xTB", xs, pick ( "xtb" ), "c1", False, True ), ( "rule", xs, pick ( "dyff_rule" ), "c2", True, False ),
                                       ( "library", xs, pick ( "dyff_torsion" ), "c5", False, False ) ],
                                     "dihedral {} (degrees)".format ( " – ".join ( str ( a + 1 ) for a in c["quad"] ) ), "relative energy (kcal/mol)",
                                     xticks = [ -180, -120, -60, 0, 60, 120, 180 ] ) )
            out.append ( '<div class="tw"><table><tr><th class="l">RMSE vs xTB rule / library</th><td class="l">{:.2f} / {:.2f} kcal/mol</td></tr>'
                         '<tr><th class="l">xTB barrier</th><td class="l">{:.2f} kcal/mol</td></tr><tr><th class="l">scan system</th><td class="l">{} atoms, charge {}, gas phase</td></tr>'
                         '{}</table></div></div>'.format ( c["rule_rmse"], c["torsion_rmse"], c["barrier_xtb"], c["fragment_atoms"], c["fragment_charge"],
                         '<tr><th class="l">logs</th><td class="l"><a href="logs/torsions/{0}/">logs/torsions/{0}/</a></td></tr>'.format ( E ( c["log_dir"] ) ) if c.get ( "log_dir" ) else "" ) )
        out.append ( "</div>" )
    details = rep.get ( "scan_details" ) or [ ]
    if "scan_targets" not in rep:
        out.append ( "<p>Automatic xTB torsion scans were not requested.</p>" ); return "\n".join ( out )
    ss = rep.get ( "scan_settings" ) or { }
    out.append ( "<h3>Scanned torsions (relaxed GFN2-xTB scan, self-consistent cosine fit)</h3>" )
    out.append ( "<p>{} flexible bond(s) without a library torsion. Each one is scanned {} in {:g}° steps; the DYFF profile is relaxed at every point "
                 "with a 1000 kJ/mol/rad² restraint. A fit is kept only when its RMSE is below {:.0%} of the DYFF rule's; otherwise the rule stays.</p>".format (
                     rep.get ( "scan_targets", 0 ),
                     "on its own capped fragment ({} bonds around it{})".format ( ss.get ( "fragment_shell", 2 ), ", gas phase" if ss.get ( "gas_phase" ) else "" ) if ss.get ( "fragments", True ) else "on the subsystem of its group",
                     ss.get ( "step", 15.0 ), ss.get ( "min_improvement", 0.9 ) ) )
    if details:
        out.append ( '<div class="tw"><table><tr><th class="l">bond</th><th class="l">dihedral</th><th>fragment</th><th>xTB barrier</th><th>RMSE rule</th><th>RMSE fit</th><th class="l">result</th></tr>' )
        for d in details:
            j, k = ( int ( x ) for x in d["bond"].split ( "-" ) )
            status = ( '<span class="ok">fitted</span>' + ( " (RMSE &gt; {} kcal/mol)".format ( ss.get ( "poor_rmse", 2.0 ) ) if d.get ( "poor" ) else "" ) ) if d.get ( "accepted" ) else \
                     ( '<span class="no">failed</span>: ' + E ( d["error"] ) if d.get ( "error" ) else '<span class="no">rejected</span>, rule kept' )
            out.append ( '<tr><td class="l"><a href="#t{}-{}">{}</a></td><td class="l">{}</td><td>{} atoms</td><td>{}</td><td>{}</td><td>{}</td><td class="l">{}</td></tr>'.format (
                j + 1, k + 1, E ( "{} – {}".format ( _alabel ( data, j ), _alabel ( data, k ) ) ), E ( " – ".join ( _alabel ( data, a ) for a in d.get ( "quad", [ ] ) ) ),
                d.get ( "fragment_atoms", "" ), _f ( d.get ( "barrier_xtb" ), "{:.2f}" ), _f ( d.get ( "before_rmse" ), "{:.2f}" ), _f ( d.get ( "fit_rmse" ), "{:.2f}" ), status ) )
        out.append ( "</table></div><p class=\"small\">Energies in kcal/mol, relative to each profile's minimum.</p>" )
    for f in rep.get ( "scan_failed" ) or [ ]:
        if not any ( f.startswith ( d["bond"] + ":" ) for d in details ): out.append ( '<p class="small">Not scanned: {}</p>'.format ( E ( f ) ) )
    out.append ( '<div class="legend"><span><i class="k1"></i>GFN2-xTB</span><span><i class="k2"></i>DYFF rule</span>'
                 '<span><i class="k3"></i>DYFF fitted (kept)</span><span><i class="k4"></i>DYFF fitted (rejected)</span></div>' )
    out.append ( '<div class="grid2">' )
    for d in details:
        j, k = ( int ( x ) for x in d["bond"].split ( "-" ) )
        p = d.get ( "profiles" ) or { }
        out.append ( '<div class="torsion" id="t{}-{}"><h3>{} – {}</h3>'.format ( j + 1, k + 1, E ( _alabel ( data, j ) ), E ( _alabel ( data, k ) ) ) )
        if p.get ( "phi" ):
            phi = [ ( ( x + 180.0 ) % 360.0 ) - 180.0 for x in p["phi"] ]
            order = sorted ( range ( len ( phi ) ), key = lambda i: phi[i] )
            xs = [ phi[i] for i in order ]
            pick = lambda key: [ p[key][i] for i in order ]
            series = [ ( "xTB", xs, pick ( "xtb" ), "c1", False, True ), ( "rule", xs, pick ( "dyff_before" ), "c2", True, False ),
                       ( "fit", xs, pick ( "dyff_after" ), "c3" if d.get ( "accepted" ) else "c4", not d.get ( "accepted" ), False ) ]
            out.append ( svg_lines ( series, "dihedral {} (degrees)".format ( " – ".join ( str ( a + 1 ) for a in d.get ( "quad", [ ] ) ) ),
                                     "relative energy (kcal/mol)", xticks = [ -180, -120, -60, 0, 60, 120, 180 ] ) )
        elif d.get ( "error" ):
            out.append ( '<p class="no">{}</p>'.format ( E ( d["error"] ) ) )
        rows = [ ( "RMSE rule / fit", "{} / {} kcal/mol".format ( _f ( d.get ( "before_rmse" ), "{:.2f}" ), _f ( d.get ( "fit_rmse" ), "{:.2f}" ) ) ),
                 ( "result", "fitted and kept" if d.get ( "accepted" ) else "not kept: DYFF rule used" ),
                 ( "scan system", "{} atoms, charge {}, {}".format ( d.get ( "fragment_atoms", "?" ), d.get ( "fragment_charge", "?" ),
                                                                     "gas phase" if not d.get ( "scan_solvent" ) else "ALPB(" + str ( d["scan_solvent"] ) + ")" ) ) ]
        if d.get ( "pairs" ):
            rows.append ( ( "fitted series (c kJ/mol × n)", "  ".join ( "{:+.3f}×{}".format ( c, int ( n ) ) for c, n in d["pairs"] ) ) )
        if d.get ( "dropped_periods" ): rows.append ( ( "periods dropped (cancel out)", ", ".join ( str ( x ) for x in d["dropped_periods"] ) ) )
        if d.get ( "log_dir" ):
            rows.append ( ( "logs", '<a href="logs/torsions/{0}/">logs/torsions/{0}/</a>'.format ( E ( d["log_dir"] ) ) ) )
        out.append ( '<div class="tw"><table>' + "".join ( '<tr><th class="l">{}</th><td class="l">{}</td></tr>'.format ( E ( a ), b if a == "logs" else E ( b ) ) for a, b in rows ) + "</table></div></div>" )
    out.append ( "</div>" )
    return "\n".join ( out )


def _logs_html ( folder ):
    E = html.escape
    root = os.path.join ( folder, "logs" )
    out = [ '<h2 id="logs">Original logs</h2>' ]
    files = [ ]
    for dirpath, dirnames, filenames in os.walk ( root ):
        dirnames.sort ( )
        for name in sorted ( filenames ): files.append ( os.path.relpath ( os.path.join ( dirpath, name ), folder ) )
    if not files:
        out.append ( "<p>No logs were kept.</p>" ); return "\n".join ( out )
    out.append ( "<p>Inputs and outputs of the external programs, as they wrote them. Click a file to open it; the outputs are also shown here.</p>" )
    current = None
    for rel in files:
        d = os.path.dirname ( rel )
        if d != current:
            if current is not None: out.append ( "</div>" )
            out.append ( '<h3>{}</h3><div>'.format ( E ( d ) ) ); current = d
        size = os.path.getsize ( os.path.join ( folder, rel ) )
        link = '<a href="{0}">{1}</a> <span class="small">({2:.1f} kB)</span>'.format ( E ( rel ), E ( os.path.basename ( rel ) ), size / 1024.0 )
        if rel.endswith ( ".out" ) and size <= CMAX_EMBED:
            with open ( os.path.join ( folder, rel ), errors = "replace" ) as handle: text = handle.read ( )
            out.append ( "<details><summary>{}</summary><pre>{}</pre></details>".format ( link, E ( text ) ) )
        else:
            out.append ( "<div>{}</div>".format ( link ) )
    out.append ( "</div>" )
    return "\n".join ( out )
