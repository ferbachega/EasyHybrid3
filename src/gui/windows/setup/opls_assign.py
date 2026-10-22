#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  opls_assign.py
#
#  [EN] 2026-10-17, user's request: assigning OPLS (System > Add > Force
#  Field > OPLS, or closing the Builder with OPLS chosen in Force Field
#  Settings) first tries the OPLS-AA parameter set as it is. Only when
#  atoms are left untyped or terms are missing, a dialog offers the
#  procedure of Extras > OPLS > Prepare OPLS System to create the new
#  parameters (pdynamo/opls/prep.prepare_opls_system: standard residues
#  from OPLS-AA, other molecules as ligands -- OPLS-AA classes by rules
#  for Lennard-Jones and bonded terms, CM5 charges from one GFN1-xTB
#  calculation, missing terms estimated; Dodda et al., JCTC 2015).
#
#  Solvation is NOT part of this: it is a later, separate step (Builder
#  Solvate, or Prepare OPLS System's own options).
#
#  When the prepared system has the same atoms in the same order (a
#  Builder molecule, a complete system) it replaces the original under the
#  same e_id, so the object stays; otherwise (hydrogens, termini or links
#  were added) it is added as a new system and the original is untouched.
#

import os
import threading
import traceback

import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk, GLib

DEFAULT_CM5_SCALE = 1.00


def default_cm5_scale ( main ):
    """ The OPLS CM5 scale of the Builder's Force Field Settings (shared). """
    settings = main.vm_session.vm_config.gl_parameters.get ( 'builder_ff_settings' ) or { }
    try:
        return float ( settings.get ( "opls_cm5_scale", DEFAULT_CM5_SCALE ) )
    except ( TypeError, ValueError ):
        return DEFAULT_CM5_SCALE


def _residue_list ( rows ):
    """ [EN] 2026-10-17, user's request: the residues that will be parametrized,
        shown in the dialog -- name, copies, atoms, and whether it is an edited
        variant of another residue (prep.split_ligand_variants()). """
    store = Gtk.ListStore ( str, str, str, str )
    for row in rows:
        note = "edited copy of {} -- parametrized separately".format ( row["original"] ) if row["original"] != row["name"] else ""
        residues = ", ".join ( row["residues"][:4] ) + ( ", ..." if len ( row["residues"] ) > 4 else "" )
        store.append ( [ row["name"], str ( row["copies"] ), "{} ({} heavy)".format ( row["atoms"], row["heavy"] ),
                         note or residues ] )
    tree = Gtk.TreeView ( model = store )
    tree.set_headers_visible ( True )
    for k, title in enumerate ( ( "Residue", "Copies", "Atoms", "Residues / note" ) ):
        renderer = Gtk.CellRendererText ( )
        column = Gtk.TreeViewColumn ( title, renderer, text = k )
        column.set_resizable ( True )
        tree.append_column ( column )
    scrolled = Gtk.ScrolledWindow ( )
    scrolled.set_policy ( Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC )
    scrolled.set_min_content_height ( min ( 220, 28 * ( len ( rows ) + 1 ) + 8 ) )
    scrolled.set_shadow_type ( Gtk.ShadowType.IN )
    scrolled.add ( tree )
    return scrolled


def _ask_parametrize ( main, system, missing, cm5_scale, parent = None, rows = None ):
    """ The question; returns the CM5 scale to use, or None (No). rows: the
        residues to parametrize (prep.ligands_to_parametrize()). """
    if missing.get ( "atoms" ):
        residues = missing.get ( "residues" ) or [ ]
        what = "{} atom(s) of '{}' have no OPLS-AA atom type (residues: {}{}).".format (
            missing["atoms"], system.label, ", ".join ( residues[:6] ), ", ..." if len ( residues ) > 6 else "" )
    else:
        what = "{} bonded term(s) of '{}' have no OPLS-AA parameters.".format ( missing.get ( "terms", 0 ), system.label )
    dialog = Gtk.MessageDialog ( transient_for = parent or main.window, modal = True,
                                 message_type = Gtk.MessageType.QUESTION, buttons = Gtk.ButtonsType.NONE,
                                 text = "Create the missing OPLS parameters?" )
    if rows:
        n_variants = sum ( 1 for r in rows if r["original"] != r["name"] )
        what = "{} residue type(s) of '{}' have no OPLS-AA parameters and will be parametrized{}:".format (
            len ( rows ), system.label, " ({} of them edited copies of another residue, renamed)".format ( n_variants ) if n_variants else "" )
    dialog.format_secondary_text (
        what + "\n\nYes: the Prepare OPLS System procedure parametrizes them -- standard residues keep the OPLS-AA "
        "parameters; the other molecules get OPLS-AA classes by rules (Lennard-Jones and bonded terms), CM5 charges "
        "from one GFN1-xTB calculation and estimated missing terms (listed in the report). Runs in the background.\n\n"
        "Solvation is not done here: solvate afterwards as a separate step.\n\n"
        "No: the system is kept without a force field." )
    if rows:
        listing = _residue_list ( rows )
        listing.show_all ( )
        dialog.get_message_area ( ).pack_start ( listing, True, True, 0 )
        dialog.set_resizable ( True )
    box = Gtk.Box ( orientation = Gtk.Orientation.HORIZONTAL, spacing = 6 )
    box.pack_start ( Gtk.Label ( label = "CM5 charge scale:" ), False, False, 0 )
    spin = Gtk.SpinButton.new_with_range ( 0.5, 1.5, 0.01 )
    spin.set_digits ( 2 )
    spin.set_value ( cm5_scale )
    spin.set_tooltip_text ( "Dodda et al., JCTC 2015 (DFT CM5 + OPLS-AA): 1.14 non-aqueous, 1.20 universal, 1.27 water. "
                            "xTB CM5 is more polar than DFT CM5: 1.00 is the EasyHybrid default. Ions keep 1.00." )
    box.pack_start ( spin, False, False, 0 )
    box.show_all ( )
    dialog.get_message_area ( ).pack_start ( box, False, False, 0 )
    dialog.add_button ( "No", Gtk.ResponseType.NO )
    dialog.add_button ( "Yes, parametrize", Gtk.ResponseType.YES )
    dialog.set_default_response ( Gtk.ResponseType.YES )
    response = dialog.run ( )
    scale = round ( spin.get_value ( ), 2 )
    dialog.destroy ( )
    return scale if response == Gtk.ResponseType.YES else None


def assign_opls ( main, system, cm5_scale = None, after = None, parent = None ):
    """ OPLS on `system`: the OPLS-AA set as it is; if parameters are
        missing, asks and runs the Extras > OPLS parametrization (no
        solvation). after( final_system ): optional, called on success (the
        same e_id system, or the new one) and may return a note for the
        message (e.g. the Builder's nonbonded model). """
    p_session = main.p_session
    if system is None:
        main.simple_dialog.error ( msg = "No active system." )
        return
    if getattr ( system, "qcModel", None ) is not None:
        main.simple_dialog.error ( msg = "'{}' has a QC model: OPLS was not assigned.".format ( system.label ) )
        return
    ok, message = p_session.define_MMModel ( force_field = 'OPLS', system = system )
    if ok:
        if after is not None:
            message += after ( system ) or ""
        main.refresh_main_statusbar ( )
        main.bottom_notebook.status_teeview_add_new_item ( message = message, system = system )
        return
    missing = getattr ( p_session, "last_opls_missing", None )
    if not missing:
        main.simple_dialog.error ( msg = message )
        return
    default = default_cm5_scale ( main ) if cm5_scale is None else cm5_scale
    main.statusbar_main.push ( 1, "OPLS: analysing the residues of '{}'...".format ( system.label ) )

    def work ( ):
        from pdynamo.opls import prep
        try:
            analysis = prep.analyze_system ( system )
            rows = prep.ligands_to_parametrize ( analysis )
            GLib.idle_add ( ask, analysis, rows, None )
        except Exception as error:
            traceback.print_exc ( )
            GLib.idle_add ( ask, None, None, str ( error ) )

    def ask ( analysis, rows, error ):
        if error is not None:
            main.simple_dialog.error ( msg = "The residues of '{}' could not be analysed:\n{}".format ( system.label, error[-600:] ) )
            return False
        scale = _ask_parametrize ( main, system, missing, default, parent, rows = rows )
        if scale is not None:
            _parametrize ( main, system, scale, after, analysis = analysis )
        else:
            main.statusbar_main.push ( 1, "OPLS: not assigned." )
        return False

    threading.Thread ( target = work, daemon = True ).start ( )


def ask_ignore_close_contact ( main, error, parent = None ):
    """ [EN] 2026-10-17, user's request: a prep.CloseContactError (two atoms of
        different residues closer than a bond length, no bond drawn between them)
        can be ignored -- only the bonds shown in the structure are bonds. True =
        ignore and continue. """
    dialog = Gtk.MessageDialog ( transient_for = parent or main.window, modal = True,
                                 message_type = Gtk.MessageType.WARNING, buttons = Gtk.ButtonsType.NONE,
                                 text = "Two atoms are very close" )
    dialog.format_secondary_text (
        "{} and {} are {:.2f} Å apart -- as close as a covalent bond -- and the procedure cannot tell whether "
        "they are bonded.\n\nIgnore and continue: only the bonds shown in the structure are bonds; these two atoms "
        "are treated as a close contact (each molecule is parametrized on its own).\n\nCancel: fix the "
        "structure first (move the atoms apart, or draw the bond if it is real).".format (
            error.atom, error.other, error.distance ) )
    dialog.add_button ( "Cancel", Gtk.ResponseType.CANCEL )
    dialog.add_button ( "Ignore and continue", Gtk.ResponseType.YES )
    dialog.set_default_response ( Gtk.ResponseType.YES )
    response = dialog.run ( )
    dialog.destroy ( )
    return response == Gtk.ResponseType.YES


def _parametrize ( main, system, cm5_scale, after, trust_connectivity = False, analysis = None ):
    p_session = main.p_session
    e_id = system.e_id
    folder = getattr ( system, "e_working_folder", None ) or os.getcwd ( )
    safe = "".join ( c if c.isalnum ( ) or c in "_-." else "_" for c in ( system.label or "system" ) )
    work_dir = os.path.join ( folder, "{}_{}_opls".format ( safe, e_id ) )
    xtb = main.vm_session.vm_config.gl_parameters.get ( 'xtb_command' ) or None
    nb_model = getattr ( system, "nbModel", None )

    def say ( text ):
        GLib.idle_add ( lambda: main.statusbar_main.push ( 1, "OPLS: " + text ) and False )
    say ( "parametrizing '{}'...".format ( system.label ) )

    def work ( ):
        from pdynamo.opls import prep
        try:
            new_system, report = prep.prepare_opls_system ( system, work_dir = work_dir, label = system.label, solvation = None,
                                                            progress = say, cm5_scale = cm5_scale, xtb = xtb, analysis = analysis,
                                                            trust_connectivity = trust_connectivity )
            text = prep.report_markdown ( report, title = "OPLS parametrization: {}".format ( system.label ) )
            GLib.idle_add ( done, new_system, report, text, None )
        except prep.CloseContactError as error:
            GLib.idle_add ( close_contact, error )
        except Exception as error:
            traceback.print_exc ( )
            GLib.idle_add ( done, None, None, None, str ( error ) )

    def close_contact ( error ):
        if ask_ignore_close_contact ( main, error ):
            _parametrize ( main, system, cm5_scale, after, trust_connectivity = True, analysis = analysis )
        else:
            main.statusbar_main.push ( 1, "OPLS: cancelled ({} / {} too close).".format ( error.atom, error.other ) )
        return False

    def done ( new_system, report, text, error ):
        failed = "The OPLS parametrization of '{}' failed".format ( system.label )
        if error is not None:
            main.simple_dialog.error ( msg = "{}:\n{}\nThe system was kept without a force field.".format ( failed, error[-600:] ) )
            return False
        report_path = os.path.join ( work_dir, "prepare_opls_report.md" )
        try:
            os.makedirs ( work_dir, exist_ok = True )
            with open ( report_path, "w" ) as handle: handle.write ( text )
        except OSError:
            report_path = None
        if not report["ok"]:
            main.simple_dialog.error_details ( msg = "{}: it could not be fully typed/parametrized.\n"
                                               "The system was kept without a force field.".format ( failed ),
                                               details = text, title = "OPLS" )
            return False
        input_files = { 'opls_parameters': report.get ( "parameter_folder" ), 'report': report_path }
        same_atoms = [ a.atomicNumber for a in new_system.atoms ] == [ a.atomicNumber for a in system.atoms ]
        if same_atoms and p_session.psystem.get ( e_id ) is system:
            # . exact original coordinates (the PDB route rounds them to 1e-3); the user's NB model kept
            for i in range ( len ( system.atoms ) ):
                for k in range ( 3 ):
                    new_system.coordinates3[i, k] = system.coordinates3[i, k]
            if nb_model is not None:
                new_system.DefineNBModel ( nb_model )
            from gui.windows.builder.empty_object import _swap_rebuilt_system_into_psystem
            _swap_rebuilt_system_into_psystem ( p_session, e_id, new_system )
            new_system.e_input_files = input_files
            message = "MM model 'OPLS-AA' assigned to '{}' with new parameters.".format ( new_system.label )
        else:
            # . hydrogens/termini/links were added: a new system, the original untouched
            name = "{}_OPLS".format ( system.label )
            try:
                new_system = p_session.add_prepared_pdynamo_system ( new_system, name = name, tag = 'OPLS', working_folder = folder,
                                                                     input_files = input_files )
            except Exception as exc:
                traceback.print_exc ( )
                main.simple_dialog.error ( msg = "The system was parametrized but could not be added:\n{}".format ( exc ) )
                return False
            message = ( "OPLS-AA with new parameters: '{}' was added as a new system ({} atoms; the preparation added "
                        "hydrogens, termini or links, so '{}' was left as it was).".format ( name, len ( new_system.atoms ), system.label ) )
        if report.get ( "ligands" ):
            message += " Parametrized: {} (CM5 x{:.2f}; {} estimated term(s)).".format (
                ", ".join ( report["ligands"] ), cm5_scale, len ( report.get ( "estimated_terms", [ ] ) ) )
        if report.get ( "ligand_variants" ):
            message += " Edited copies renamed: {}.".format ( ", ".join ( "{} (from {})".format ( new, old )
                                                                   for new, old in report["ligand_variants"].items ( ) ) )
        for path, data in ( report.get ( "quality" ) or { } ).items ( ):
            levels = data["levels"]
            message += " Quality {}: {} green, {} yellow, {} red.".format ( path, levels["green"], levels["yellow"], levels["red"] )
        if after is not None:
            message += after ( new_system ) or ""
        message += " Report: {}".format ( report_path ) if report_path else ""
        main.refresh_main_statusbar ( )
        main.main_treeview.refresh ( )
        main.bottom_notebook.status_teeview_add_new_item ( message = message, system = new_system )
        main.statusbar_main.push ( 1, message )
        # . the quality report of the new parameters (Force Field Parameters window)
        window = getattr ( main, "opls_parameters_window", None )
        if window is not None:
            try:
                window.open_window ( e_id = getattr ( new_system, "e_id", None ) )
            except Exception:
                traceback.print_exc ( )
        return False

    threading.Thread ( target = work, daemon = True ).start ( )
