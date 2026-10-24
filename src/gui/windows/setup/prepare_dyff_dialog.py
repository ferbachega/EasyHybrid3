#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: "Prepare DYFF System" dialog
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description (2026-10-08):
#      A raw structure (PDB without bonds) -> pdynamo/dyff_prep.py (hydrogens,
#      protonation, metal sites, truncated residues, ligands; no OPLS typing) -> a new
#      system in the session -> DYFF + residue library applied with the options chosen
#      here (main_window._run_dyff_library). Opened by System > Add > Force Field >
#      DYFF + residue library... when the active system has no bonds, and by
#      Extras > Prepare DYFF System...
#

import os
import threading
import traceback

import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import GLib, Gtk


#  2026-10-21, user's request: System > Force Field > Define Force Field > Prepare
#  DYFF System first asks which DYFF to use.
DYFF_VERSIONS = (
    ( "modified", "DYFF modified: residue library + Seminario + CM5 (recommended)",
      "Residues matched to the precomputed xTB library; the rest refined on capped fragments: bonds and angles "
      "from the xTB Hessian (modified Seminario), torsions from the DFT torsion library, CM5 charges. A raw "
      "structure (PDB without bond orders) is prepared first: hydrogens, protonation, metal sites, ligands." ),
    ( "eh", "DYFF 1.0-EH: generic rules, EasyHybrid atom types",
      "pDynamo's DYFF rules (UFF-derived) with the atom types EasyHybrid adds (aromatic P/Se/B, boranes, "
      "metals). Needs bond orders. No charges are computed: partial charges already on the object (Builder, "
      "a previous residue-library run) are used, otherwise all charges are zero." ),
    ( "pdynamo", "DYFF 1.0: pDynamo's original parameter set",
      "The dyff-1.0 set shipped with pDynamo3, unchanged. Same rules and charge handling as above; atoms "
      "outside its types (aromatic P/Se, boron clusters...) are not typed." ),
)


def choose_dyff_version ( main ):
    """ ( DYFF version chosen ("modified", "eh", "pdynamo"), new system? ) or None.
        2026-10-21, user's request: a new system (a copy) is the default. """
    dialog = Gtk.Dialog ( title = "Prepare DYFF System", transient_for = main.window, modal = True )
    dialog.add_buttons ( "Cancel", Gtk.ResponseType.CANCEL, "Continue", Gtk.ResponseType.OK )
    dialog.set_default_response ( Gtk.ResponseType.OK )
    box = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 10, margin = 12 )
    head = Gtk.Label ( xalign = 0 ); head.set_markup ( "<b>Which DYFF version?</b>" )
    box.pack_start ( head, False, False, 0 )
    group, buttons = None, { }
    for key, title, text in DYFF_VERSIONS:
        button = Gtk.RadioButton.new_with_label_from_widget ( group, title )
        group = group or button
        detail = Gtk.Label ( label = text, xalign = 0, wrap = True, margin_start = 26 )
        detail.set_max_width_chars ( 80 ); detail.get_style_context ( ).add_class ( "dim-label" )
        box.pack_start ( button, False, False, 0 ); box.pack_start ( detail, False, False, 0 )
        buttons[key] = button
    box.pack_start ( Gtk.Separator ( ), False, False, 2 )
    head = Gtk.Label ( xalign = 0 ); head.set_markup ( "<b>Result</b>" )
    box.pack_start ( head, False, False, 0 )
    radio_new  = Gtk.RadioButton.new_with_label_from_widget ( None, "New system (a copy named <system>_DYFF; the original is kept)" )
    radio_here = Gtk.RadioButton.new_with_label_from_widget ( radio_new, "Apply to this system" )
    note = Gtk.Label ( label = "A raw structure (PDB without bonds) is rebuilt by \"DYFF modified\" and always gives a new system.",
                       xalign = 0, wrap = True, margin_start = 26 )
    note.set_max_width_chars ( 80 ); note.get_style_context ( ).add_class ( "dim-label" )
    for w in ( radio_new, radio_here, note ): box.pack_start ( w, False, False, 0 )
    dialog.get_content_area ( ).add ( box )
    dialog.show_all ( )
    response = dialog.run ( )
    chosen = next ( ( k for k, b in buttons.items ( ) if b.get_active ( ) ), None )
    as_new = radio_new.get_active ( )
    dialog.destroy ( )
    return ( chosen, as_new ) if response == Gtk.ResponseType.OK else None


def prepare_dyff_system ( main ):
    """ Menu entry: version choice, then the matching path. """
    system = main.p_session.psystem.get ( main.p_session.active_id )
    if system is None:
        main.simple_dialog.error ( msg = "No active system." ); return
    if getattr ( system, "qcModel", None ) is not None:
        main.simple_dialog.error ( msg = "The system has a QC model: remove the QC region before redefining the MM model." ); return
    choice = choose_dyff_version ( main )
    if choice is None: return
    version, as_new = choice
    if version == "modified":
        main._dyff_residue_library_dialog ( as_new = as_new ); return
    from pdynamo.dyff_edits import EH_DYFF_SET, PDYNAMO_DYFF_SET
    if as_new:
        vobject = next ( ( v for v in main.vm_session.vm_objects_dic.values ( ) if getattr ( v, "e_id", None ) == system.e_id ), None )
        system, _ = main._dyff_copy_system ( system, vobject )
    # . the per-term edits of a previous residue-library run do not belong to the generic rules
    for vobject in main.vm_session.vm_objects_dic.values ( ):
        if getattr ( vobject, "e_id", None ) == system.e_id:
            vobject.dyff_parameter_edits = None
            vobject.dyff_library_options = None
    ok, msg = main.p_session.define_MMModel ( force_field = "DYFF", system = system,
                                              parameter_set = EH_DYFF_SET if version == "eh" else PDYNAMO_DYFF_SET )
    main.refresh_main_statusbar ( )
    if ok:
        main.bottom_notebook.status_teeview_add_new_item ( message = msg, system = system )
    else:
        main.simple_dialog.error ( msg = msg + "\n\nA structure without bond orders (raw PDB) needs \"DYFF modified\", which prepares it first." )


def _headless_charge_multiplicity ( system, options ):
    """ ( charge, multiplicity or None ) for the xTB calculations without asking: the
        defaults of the charge/multiplicity dialog (sum of the formal charges; the metal
        sites keep their own multiplicity) unless options give "charge" / "mult". A
        multiplicity that does not agree with the electron count -> error text. """
    electrons = sum ( int ( a.atomicNumber ) for a in system.atoms )
    q = options.get ( "charge" )
    q = sum ( int ( getattr ( a, "formalCharge", 0 ) or 0 ) for a in system.atoms ) if q is None else int ( q )
    m = options.get ( "mult" )
    if m is None:
        return q, None
    m = int ( m )
    if ( electrons - q ) % 2 != ( m - 1 ) % 2:
        return ( "charge {:+d} gives {} electrons: the multiplicity must be {}.".format (
                 q, electrons - q, "odd" if ( electrons - q ) % 2 == 0 else "even" ) )
    return q, m


def prepare_dyff_headless ( main, system, version = "modified", new = True, name = None, ph = 7.4, water = True,
                            method = "gfn2", scale = 1.0, scan = False, complete = True, remove_caps = True,
                            charge = None, mult = None, ligand_charges = "cm5" ):
    """ 2026-10-22: Prepare DYFF System without dialogs (terminal command 'dyff').
        ligand_charges (2026-10-23): cm5 | resp-xtb | resp-hf (residue_library.LIGAND_CHARGE_MODELS).
        Same paths as the menu entry: "eh"/"pdynamo" define DYFF from the bond orders;
        "modified" = residue library + xTB for the rest -- a raw structure (no bond
        orders) is built first by dyff_prep (hydrogens, protonation, metal sites with
        their default charge/multiplicity, caps removed, truncated residues completed)
        and always gives a new system. charge/mult: of the xTB calculations (default:
        sum of the formal charges; the metal sites' multiplicity). Runs in the background;
        the result is reported in the Status panel. Returns a text for the terminal. """
    if system is None:
        return "No system."
    if getattr ( system, "qcModel", None ) is not None:
        return "The system has a QC model: remove the QC region before redefining the MM model."
    if method not in ( "gfn2", "gxtb" ):
        return "method must be gfn2 or gxtb."
    if name and any ( s is not None and s.label == name for s in main.p_session.psystem.values ( ) ):
        return "A system named '{}' already exists: choose another name=.".format ( name )
    vobject = next ( ( v for v in main.vm_session.vm_objects_dic.values ( ) if getattr ( v, "e_id", None ) == system.e_id ), None )
    if version in ( "eh", "pdynamo" ):
        from pdynamo.dyff_edits import EH_DYFF_SET, PDYNAMO_DYFF_SET
        if new:
            system, vobject = main._dyff_copy_system ( system, vobject )
            if name:
                main.rename ( e_id = system.e_id, name = name )
        if vobject is not None:
            vobject.dyff_parameter_edits = None
            vobject.dyff_library_options = None
        ok, msg = main.p_session.define_MMModel ( force_field = "DYFF", system = system,
                                                  parameter_set = EH_DYFF_SET if version == "eh" else PDYNAMO_DYFF_SET )
        main.refresh_main_statusbar ( )
        if ok:
            main.bottom_notebook.status_teeview_add_new_item ( message = msg, system = system )
            return msg
        return msg + " A structure without bond orders (raw PDB) needs version=modified."
    if version != "modified":
        return "version must be modified, eh or pdynamo."
    options = { "complete": bool ( complete ), "ph": float ( ph ), "water": bool ( water ), "method": method,
                "scale": float ( scale ), "scan_torsions": bool ( scan ), "headless": True, "charge": charge, "mult": mult,
                "ligand_charges": ligand_charges,
                "name": name or "{}_DYFF".format ( system.label or "system" ) }
    if main._dyff_needs_preparation ( system ):
        from pdynamo import dyff_prep as DP
        info = DP.analyze ( system )
        if remove_caps and info["covalent"]:
            DP.remove_residues ( info, list ( info["covalent"] ) )
        PrepareDYFFDialog ( main )._build ( system, info, options )
        sites = ", ".join ( "{} {}{:+d} {} (mult {})".format ( s["label"], s["element"], s["charge"], s["model"], s["multiplicity"] )
                            for s in info["metals"] ) or "none"
        return ( "Prepare DYFF System started: new system '{}' (raw structure: hydrogens, protonation at pH {:.1f}, metal sites: {}). "
                 "xTB runs in the background; the result appears in the Status panel.".format ( options["name"], options["ph"], sites ) )
    # . bond orders already known: residue library directly
    from pdynamo import residue_library as RL
    if vobject is None or len ( vobject.atoms ) != len ( system.atoms ):
        return "No object with the same atoms as the system."
    if new:
        system, vobject = main._dyff_copy_system ( system, vobject )
        if vobject is None:
            return "The copy of the system has no object."
    if new and name:
        main.rename ( e_id = system.e_id, name = name )
    p = RL.plan ( system )
    q, m = None, None
    if p["summary"]["rest_atoms"]:
        answer = _headless_charge_multiplicity ( system, options )
        if isinstance ( answer, str ):
            return answer
        q, m = answer
    main._run_dyff_library ( system, vobject, p, use_water = options["water"], use_bonded = True, method = method,
                             cm5_scale = options["scale"], scan_torsions = options["scan_torsions"], total_charge = q, multiplicity = m,
                             ligand_charges = ligand_charges )
    s = p["summary"]
    return ( "DYFF modified started on '{}': {}/{} residues from the library, {} atoms with xTB in the background; "
             "the result appears in the Status panel.".format ( system.label, s["matched"], s["residues"], s["rest_atoms"] ) )


class PrepareDYFFDialog:

    def __init__ ( self, main ):
        self.main = main

    def run ( self, system ):
        from pdynamo import dyff_prep as DP
        try:
            info = DP.analyze ( system )
        except Exception as error:
            traceback.print_exc ( )
            self.main.simple_dialog.error ( msg = "Prepare DYFF System: the structure could not be analyzed:\n{}".format ( error ) ); return
        analysis = info["analysis"]
        residues = analysis.residues
        kinds = { }
        for r in residues: kinds[r.kind] = kinds.get ( r.kind, 0 ) + 1

        dialog = Gtk.Dialog ( title = "Prepare DYFF System", transient_for = self.main.window, modal = True )
        dialog.add_buttons ( "Cancel", Gtk.ResponseType.CANCEL, "Prepare and apply DYFF", Gtk.ResponseType.OK )
        dialog.set_default_size ( 820, 640 )
        box = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 8, margin = 10 )
        scroll = Gtk.ScrolledWindow ( ); scroll.set_vexpand ( True ); scroll.add ( box )
        dialog.get_content_area ( ).pack_start ( scroll, True, True, 0 )
        def label ( text, bold = False, wrap = True ):
            w = Gtk.Label ( xalign = 0 ); w.set_line_wrap ( wrap ); w.set_max_width_chars ( 110 )
            w.set_markup ( "<b>{}</b>".format ( GLib.markup_escape_text ( text ) ) if bold else GLib.markup_escape_text ( text ) )
            return w

        box.pack_start ( label ( "{}: {} residues ({}); {} chain break(s); {} disulfide(s). No OPLS typing: hydrogens, protonation, bond "
                                 "orders, metal sites and ligands are prepared for DYFF + residue library.".format (
                                 system.label, len ( residues ), ", ".join ( "{} {}".format ( n, k ) for k, n in sorted ( kinds.items ( ) ) ),
                                 len ( analysis.chain_breaks ), len ( analysis.disulfides ) ) ), False, False, 0 )

        # . covalently bound residues the preparation cannot link (peptide caps, covalent ligands)
        remove_cov = None
        if info["covalent"]:
            box.pack_start ( label ( "Covalently bound residues without a library component", True ), False, False, 0 )
            box.pack_start ( label ( ", ".join ( residues[i].label for i in info["covalent"] ) +
                                     ": terminal caps (ACE, NME, NH2) or covalent ligands. The preparation cannot link them." ), False, False, 0 )
            remove_cov = Gtk.CheckButton ( label = "Remove them (a removed cap leaves a free N- or C-terminus)" )
            remove_cov.set_active ( True ); box.pack_start ( remove_cov, False, False, 0 )

        # . truncated residues
        complete = None
        if info["truncated"]:
            box.pack_start ( label ( "Truncated residues", True ), False, False, 0 )
            box.pack_start ( label ( ", ".join ( "{} (missing {})".format ( l, " ".join ( m ) ) for _, l, m in info["truncated"][:20] ) +
                                     ( " ..." if len ( info["truncated"] ) > 20 else "" ) ), False, False, 0 )
            complete = Gtk.CheckButton ( label = "Complete the truncated residues (WARNING: the missing atoms are placed geometrically, "
                                                 "extended conformation -- minimize and check them before production)" )
            complete.set_active ( True ); box.pack_start ( complete, False, False, 0 )

        # . metal sites (editable)
        store = None
        if info["metals"]:
            box.pack_start ( label ( "Metal sites", True ), False, False, 0 )
            box.pack_start ( label ( "Bonded: metal-donor bonds, the metal and its donor residues are refined together with xTB (they are "
                                     "\"pending\" until then: no library parameters). Ionic: no bonds, the metal is a point charge with "
                                     "Lennard-Jones (usual for Mg and Ca). Check the charge (oxidation state) and the spin multiplicity." ),
                             False, False, 0 )
            store = Gtk.ListStore ( str, str, str, int, int, int )        # site, element, donors, model index, charge, multiplicity
            for k, site in enumerate ( info["metals"] ):
                donors = ", ".join ( "{} {} {:.2f}".format ( residues[j].label.split ( ":" )[-1], n, d ) for j, n, d in site["donors"] ) or "none found"
                store.append ( [ site["label"], site["element"], donors, 0 if site["model"] == "bonded" else 1, site["charge"], site["multiplicity"] ] )
            models = Gtk.ListStore ( str ); models.append ( [ "bonded" ] ); models.append ( [ "ionic" ] )
            view = Gtk.TreeView ( model = store )
            for col, title in ( ( 0, "Site" ), ( 1, "Metal" ), ( 2, "Donors (A)" ) ):
                view.append_column ( Gtk.TreeViewColumn ( title, Gtk.CellRendererText ( ), text = col ) )
            model_cell = Gtk.CellRendererCombo ( model = models, text_column = 0, has_entry = False, editable = True )
            def model_text ( column, cell, tree_model, it, data ):
                cell.set_property ( "text", [ "bonded", "ionic" ][tree_model[it][3]] )
            def model_edited ( cell, path, text ): store[path][3] = 0 if text == "bonded" else 1
            model_cell.connect ( "edited", model_edited )
            col = Gtk.TreeViewColumn ( "Model", model_cell ); col.set_cell_data_func ( model_cell, model_text ); view.append_column ( col )
            for idx, title, lo, hi in ( ( 4, "Charge", -4, 8 ), ( 5, "Multiplicity", 1, 8 ) ):
                cell = Gtk.CellRendererSpin ( editable = True, adjustment = Gtk.Adjustment ( 0, lo, hi, 1, 1, 0 ) )
                def edited ( cell, path, text, idx = idx ):
                    try: store[path][idx] = int ( float ( text ) )
                    except ValueError: pass
                cell.connect ( "edited", edited )
                view.append_column ( Gtk.TreeViewColumn ( title, cell, text = idx ) )
            box.pack_start ( view, False, False, 0 )
        if info["metal_cofactors"]:
            names = sorted ( { residues[i].name for i in info["metal_cofactors"] } )
            box.pack_start ( label ( "Not supported yet: metal-containing cofactors {} (heme, Fe-S clusters...). Remove them or use "
                                     "Prepare OPLS System / an active-site model in the Builder.".format ( ", ".join ( names ) ) ), False, False, 0 )

        free_ligands = [ i for i in info["ligands"] if i not in set ( info["covalent"] ) ]
        if free_ligands:
            names = sorted ( { residues[i].name for i in free_ligands } )
            box.pack_start ( label ( "Ligands (hydrogens and bond orders from Open Babel at the pH below; refined like the metal sites): {}".format (
                                     ", ".join ( names ) ) ), False, False, 0 )

        grid = Gtk.Grid ( column_spacing = 8, row_spacing = 6 )
        grid.attach ( Gtk.Label ( label = "pH (ligand hydrogens):", xalign = 0 ), 0, 0, 1, 1 )
        ph = Gtk.SpinButton.new_with_range ( 0.0, 14.0, 0.1 ); ph.set_digits ( 1 ); ph.set_value ( 7.4 ); grid.attach ( ph, 1, 0, 1, 1 )
        water = Gtk.CheckButton ( label = "Implicit water (ALPB) for the CM5 charges and GFN2-xTB Hessians" ); water.set_active ( True )
        grid.attach ( water, 0, 1, 2, 1 )
        method_combo, scale_spin = self.main._dyff_method_and_scale_widgets ( grid, 2 )
        scan = self.main._dyff_scan_torsions_widget ( )
        grid.attach ( scan, 0, 5, 2, 1 )
        grid.attach ( Gtk.Label ( label = "Ligand charges:", xalign = 0 ), 0, 6, 1, 1 )
        ligand_combo = self.main._dyff_ligand_charges_widget ( ); grid.attach ( ligand_combo, 1, 6, 1, 1 )
        grid.attach ( Gtk.Label ( label = "New system name:", xalign = 0 ), 0, 4, 1, 1 )
        name_entry = Gtk.Entry ( ); name_entry.set_text ( "{}_DYFF".format ( system.label or "system" ) ); grid.attach ( name_entry, 1, 4, 1, 1 )
        box.pack_start ( grid, False, False, 0 )
        dialog.show_all ( )
        response = dialog.run ( )
        if remove_cov is not None and remove_cov.get_active ( ):
            DP.remove_residues ( info, list ( info["covalent"] ) )
        options = { "complete": complete.get_active ( ) if complete is not None else True, "ph": ph.get_value ( ),
                    "water": water.get_active ( ), "method": method_combo.get_active_id ( ), "scale": scale_spin.get_value ( ),
                    "scan_torsions": scan.get_active ( ),
                    "ligand_charges": ligand_combo.get_active_id ( ) or "cm5",
                    "name": name_entry.get_text ( ).strip ( ) or "{}_DYFF".format ( system.label ) }
        if store is not None:
            for k, row in enumerate ( store ):
                info["metals"][k]["model"] = "bonded" if row[3] == 0 else "ionic"
                info["metals"][k]["charge"], info["metals"][k]["multiplicity"] = int ( row[4] ), int ( row[5] )
        dialog.destroy ( )
        if response != Gtk.ResponseType.OK: return
        self._build ( system, info, options )

    def _build ( self, system, info, options ):
        from pdynamo import dyff_prep as DP
        self.main.statusbar_main.push ( 1, "Prepare DYFF System: building {}...".format ( options["name"] ) )
        def work ( ):
            try:
                changed = DP.apply_site_protonation ( info )
                new_system, report = DP.build ( info, label = options["name"], complete_truncated = options["complete"], ph = options["ph"] )
                report["protonation_changes"] = changed
                GLib.idle_add ( self._built, system, new_system, report, options, None )
            except Exception as error:
                traceback.print_exc ( )
                GLib.idle_add ( self._built, system, None, None, options, str ( error ) )
        threading.Thread ( target = work, daemon = True ).start ( )

    def _built ( self, source, system, report, options, error ):
        from pdynamo import dyff_prep as DP
        from pdynamo import residue_library as RL
        if error:
            self.main.simple_dialog.error ( msg = "Prepare DYFF System failed:\n{}".format ( error[-800:] ) ); return False
        if report.get ( "undefined" ):
            # . atoms left at pDynamo's undefined-coordinate sentinel (1e4 A) would end up in the
            #   object: its mass centre far from the molecule, so the view centres on empty space
            paths = report["undefined"]
            self.main.simple_dialog.error ( msg = "Prepare DYFF System: {} atoms could not be given coordinates "
                                                  "(their names do not match the residue templates), e.g.:\n{}\n\nNo system was created.".format (
                                                  len ( paths ), "\n".join ( paths[:8] ) ) )
            return False
        folder = getattr ( source, "e_working_folder", None ) or os.getcwd ( )
        try:
            system = self.main.p_session.add_prepared_pdynamo_system ( system, name = options["name"], tag = "DYFF", working_folder = folder )
        except Exception as exc:
            traceback.print_exc ( )
            self.main.simple_dialog.error ( msg = "The system was prepared but could not be added:\n{}".format ( exc ) ); return False
        vobject = next ( ( v for v in self.main.vm_session.vm_objects_dic.values ( ) if getattr ( v, "e_id", None ) == system.e_id ), None )
        if vobject is not None and report.get ( "type_overrides" ):
            # . DYFF types of the metals (dyff_prep.type_overrides), used by define_MMModel
            vobject.dyff_type_overrides = dict ( report["type_overrides"] )
        p = RL.plan ( system )
        notes = [ ]
        if report.get ( "completed" ): notes.append ( "{} atoms of truncated residues placed geometrically".format ( len ( report["completed"] ) ) )
        if report.get ( "protonation_changes" ): notes.append ( "protonation set by the metal sites: " +
                                                               ", ".join ( "{} {}".format ( l.split ( ":" )[-1], v ) for l, v in report["protonation_changes"] ) )
        bonded = [ s for s in report.get ( "metal_sites", [ ] ) if s["model"] == "bonded" ]
        if bonded: notes.append ( "{} bonded metal site(s) refined with {}".format ( len ( bonded ), { "gxtb": "g-xTB" }.get ( options["method"], "GFN2-xTB" ) ) )
        extra = ( " Prepared: " + "; ".join ( notes ) + "." ) if notes else ""
        for w in report.get ( "warnings", [ ] ): extra += " WARNING: " + w
        site_states = DP.site_groups ( system, report )
        total_charge, multiplicity = None, None
        if p["summary"]["rest_atoms"] and vobject is not None:
            # . 2026-10-21, user's request: the charge (and spin) of the parametrization; asked
            #   here because the electrons are only known once the hydrogens are built
            vobject.dyff_library_options = { "multiplicities": { system.atoms[i].path: st["multiplicity"] for i, st in site_states.items ( ) } }
            ligands = ", ".join ( "{} {:+d}".format ( l["name"], int ( l["charge"] ) ) for l in report.get ( "ligands", [ ] ) )
            intro = ( "\"{}\" was built: {} atoms, sum of the formal charges {:+d}{}. {} atoms are outside the residue library and are "
                      "parametrized with xTB (CM5 charges{}). Check the total charge and the multiplicity of these calculations.".format (
                      system.label, len ( system.atoms ), int ( report.get ( "formal_charge", 0 ) ),
                      " (ligands: {})".format ( ligands ) if ligands else "", p["summary"]["rest_atoms"],
                      ", Hessians, torsions" ) )
            if options.get ( "headless" ):
                # . terminal (prepare_dyff_headless): no dialog, the charge/multiplicity given there
                answer = _headless_charge_multiplicity ( system, options )
                if isinstance ( answer, str ):
                    self.main.simple_dialog.error ( msg = "Prepare DYFF System: " + answer ); return False
            else:
                answer = self.main._ask_dyff_charge_multiplicity ( system, vobject, intro )
            if answer is None:
                self.main.bottom_notebook.status_teeview_add_new_item (
                    message = "Prepare DYFF System: {} was built without a force field (parametrization cancelled).".format ( system.label ), system = system )
                return False
            total_charge, multiplicity = answer
            if multiplicity is not None and site_states:
                site_states = None          # . _run_dyff_library puts it on the first metal site
        self.main._run_dyff_library ( system, vobject, p, use_water = options["water"], use_bonded = True, method = options["method"],
                                      site_states = site_states, cm5_scale = options["scale"], extra_message = extra,
                                      scan_torsions = options.get ( "scan_torsions", False ),
                                      ligand_charges = options.get ( "ligand_charges", "cm5" ),
                                      total_charge = total_charge, multiplicity = multiplicity, open_report = True )
        return False
