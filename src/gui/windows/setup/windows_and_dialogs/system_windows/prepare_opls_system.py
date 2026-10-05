#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Prepare OPLS System window
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      Extras > OPLS > Prepare OPLS System. GUI over pdynamo/opls:
#        Analyze  -> prep.analyze_system() on the chosen system/frame, shown
#                    as a residue table (kind, status, editable protonation);
#        Prepare  -> prep.prepare_opls_system() in a worker thread
#                    (hydrogens, termini, disulfides, ligand parametrization
#                    with CM5 charges, solvation, OPLS typing, per-system
#                    parameter set) and the result added as a NEW system.
#      A Markdown report is written into the working folder.
#

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

import os
import threading
import traceback

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import FolderChooserButton


_VARIANT_TEXT = { None                 : "default",
                  "Delta Protonated"   : "HID (delta)",
                  "Epsilon Protonated" : "HIE (epsilon)",
                  "Doubly Protonated"  : "HIP (+1)",
                  "Deprotonated"       : "deprotonated",
                  "Protonated"         : "protonated",
                  "Neutral"            : "neutral" }

_FILTERS = [ ( "all", "All residues" ), ( "attention", "Needs attention" ),
             ( "titratable", "Titratable (HIS/ASP/GLU/LYS/CYS)" ), ( "ligands", "Ligands / non-standard" ) ]

( _COL_INDEX, _COL_LABEL, _COL_KIND, _COL_STATUS, _COL_VARIANT, _COL_EDITABLE, _COL_ATTENTION ) = range ( 7 )


class PrepareOPLSSystemWindow:
    """ "Prepare OPLS System" window. """

    def __init__ ( self, main = None ):
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        self.home       = main.home
        self.Visible    = False
        self.analysis   = None
        self.running    = False

    #-------------------------------------------------------------------------
    #  W I N D O W   L I F E C Y C L E
    #-------------------------------------------------------------------------
    def open_window ( self ):
        if self.Visible:
            self.window.present ( )
            return
        self.builder = Gtk.Builder ( )
        self.builder.add_from_file ( os.path.join ( self.home, 'src/gui/windows/setup/prepare_opls_system_window.glade' ) )
        self.builder.connect_signals ( self )
        self.window = self.builder.get_object ( 'window' )

        get = self.builder.get_object
        self.treeview          = get ( 'treeview_residues' )
        self.combo_filter      = get ( 'combo_filter' )
        self.label_count       = get ( 'label_count' )
        self.label_summary     = get ( 'label_summary' )
        self.label_status      = get ( 'label_status' )
        self.spinner           = get ( 'spinner' )
        self.button_run        = get ( 'button_run' )
        self.button_analyze    = get ( 'button_analyze' )
        self.checkbox_truncate = get ( 'checkbox_truncate' )
        self.checkbox_parametrize = get ( 'checkbox_parametrize' )
        self.spin_ph           = get ( 'spin_ph' )
        self.spin_cm5_scale    = get ( 'spin_cm5_scale' )
        self.entry_xtb         = get ( 'entry_xtb' )
        self.checkbox_solvate  = get ( 'checkbox_solvate' )
        self.grid_solvation    = get ( 'grid_solvation' )
        self.combo_shape       = get ( 'combo_shape' )
        self.spin_padding      = get ( 'spin_padding' )
        self.checkbox_align    = get ( 'checkbox_align' )
        self.checkbox_neutralize = get ( 'checkbox_neutralize' )
        self.spin_concentration = get ( 'spin_concentration' )
        self.combo_cation      = get ( 'combo_cation' )
        self.combo_anion       = get ( 'combo_anion' )
        self.spin_closeness    = get ( 'spin_closeness' )
        self.entry_system_name = get ( 'entry_system_name' )
        self.checkbox_open_parameters = get ( 'checkbox_open_parameters' )

        # . system / coordinates / folder widgets
        self.combobox_systems = SystemComboBox ( self.main )
        self.combobox_systems.connect ( "changed", self.on_combobox_systems_changed )
        get ( 'box_system' ).pack_start ( self.combobox_systems, False, False, 0 )
        self.coordinates_combobox = CoordinatesComboBox ( )
        get ( 'box_coordinates' ).pack_start ( self.coordinates_combobox, False, False, 0 )
        self.work_folder_chooser = FolderChooserButton ( self.main, 'folder', self.home )
        get ( 'box_work_folder' ).pack_start ( self.work_folder_chooser.btn, False, False, 0 )

        # . residue and ligand tables
        self._build_residue_table ( )
        self._build_ligand_table ( )
        self.spin_ph.connect ( "value-changed", lambda widget: self._fill_ligand_table ( ) )
        for key, text in _FILTERS:
            self.combo_filter.append ( key, text )
        self.combo_filter.set_active ( 0 )

        # . xtb
        try:
            from pdynamo.opls import ligand
            xtb = self.vm_session.vm_config.gl_parameters.get ( 'xtb_command' ) or ligand.find_xtb ( ) or ''
        except Exception:
            xtb = ''
        self.entry_xtb.set_text ( xtb )

        if self.p_session.psystem:
            self.combobox_systems.set_active_system ( e_id = self.p_session.active_id )
        self.window.show_all ( )
        self.spinner.hide ( )
        self.window.connect ( 'destroy', self.close_window )
        self.Visible = True

    def close_window ( self, button = None, data = None ):
        if not self.Visible: return
        self.window.destroy ( )
        self.Visible  = False
        self.analysis = None

    #-------------------------------------------------------------------------
    #  R E S I D U E   T A B L E
    #-------------------------------------------------------------------------
    def _build_residue_table ( self ):
        self.liststore = Gtk.ListStore ( int, str, str, str, str, bool, bool )
        self.filter_model = self.liststore.filter_new ( )
        self.filter_model.set_visible_func ( self._row_visible )
        self.treeview.set_model ( self.filter_model )
        for title, column in ( ( "Residue", _COL_LABEL ), ( "Kind", _COL_KIND ), ( "Status", _COL_STATUS ) ):
            renderer = Gtk.CellRendererText ( )
            col = Gtk.TreeViewColumn ( title, renderer, text = column )
            col.set_resizable ( True )
            col.set_sort_column_id ( column )
            self.treeview.append_column ( col )
        renderer = Gtk.CellRendererCombo ( )
        renderer.set_property ( "text-column", 0 )
        renderer.set_property ( "has-entry", False )
        renderer.connect ( "edited", self.on_variant_edited )
        col = Gtk.TreeViewColumn ( "Protonation", renderer, text = _COL_VARIANT, editable = _COL_EDITABLE )
        col.set_cell_data_func ( renderer, self._variant_cell_data )
        self.treeview.append_column ( col )
        self._variant_models = { }

    #-------------------------------------------------------------------------
    #  L I G A N D   T A B L E   (charge / multiplicity)
    #-------------------------------------------------------------------------
    # . columns: name, copies, atoms text, perceived charge text, charge, multiplicity,
    #            electrons with the perceived charge (-1 = unknown), perceived charge, status
    def _build_ligand_table ( self ):
        self.ligand_store = Gtk.ListStore ( str, int, str, str, int, int, int, int, str )
        self.treeview_ligands = self.builder.get_object ( 'treeview_ligands' )
        self.treeview_ligands.set_model ( self.ligand_store )
        for title, column in ( ( "Ligand", 0 ), ( "Copies", 1 ), ( "Atoms", 2 ), ( "Perceived charge", 3 ) ):
            self.treeview_ligands.append_column ( Gtk.TreeViewColumn ( title, Gtk.CellRendererText ( ), text = column ) )
        for title, column, lower, upper in ( ( "Total charge", 4, -20, 20 ), ( "Multiplicity", 5, 1, 9 ) ):
            renderer = Gtk.CellRendererSpin ( )
            renderer.set_property ( "editable", True )
            renderer.set_property ( "adjustment", Gtk.Adjustment ( value = 0, lower = lower, upper = upper, step_increment = 1, page_increment = 1 ) )
            renderer.connect ( "edited", self.on_ligand_value_edited, column )
            self.treeview_ligands.append_column ( Gtk.TreeViewColumn ( title, renderer, text = column ) )
        self.treeview_ligands.append_column ( Gtk.TreeViewColumn ( "Status", Gtk.CellRendererText ( ), text = 8 ) )

    def _ligand_status ( self, row ):
        electrons_perceived, perceived = row[6], row[7]
        if electrons_perceived < 0: return "chemistry not perceived"
        electrons = electrons_perceived - ( row[4] - perceived )
        if ( electrons + row[5] - 1 ) % 2 != 0:
            return "{} electrons: needs an {} multiplicity".format ( electrons, "odd" if electrons % 2 == 0 else "even" )
        text = "{} electrons".format ( electrons )
        if row[4] != perceived: text += "; charge changed by the user"
        return text

    def _fill_ligand_table ( self ):
        from pdynamo.opls import prep
        if not hasattr ( self, "ligand_store" ): return
        self.ligand_store.clear ( )
        if self.analysis is None or not self.analysis.unsupported: return
        try:
            preview = prep.ligand_preview ( self.analysis, ph = self.spin_ph.get_value ( ) )
        except Exception:
            traceback.print_exc ( )
            return
        for name, entry in preview.items ( ):
            perceived = entry["charge"] if entry["charge"] is not None else 0
            atoms = "{} heavy + {} H".format ( entry["heavy"], entry["hydrogens"] ) if entry["hydrogens"] is not None else "{} heavy".format ( entry["heavy"] )
            row = [ name, entry["copies"], atoms, "{:+d}".format ( perceived ) if entry["charge"] is not None else "?",
                    perceived, entry["multiplicity"], entry["electrons"] if entry["electrons"] is not None else -1, perceived, "" ]
            row[8] = entry["error"][:80] if entry["error"] else self._ligand_status ( row )
            self.ligand_store.append ( row )

    def on_ligand_value_edited ( self, renderer, path, text, column ):
        try:
            value = int ( float ( text ) )
        except ValueError:
            return
        row = self.ligand_store[path]
        if column == 5 and value < 1: value = 1
        row[column] = value
        row[8] = self._ligand_status ( row )

    def _ligand_options ( self ):
        """ { name: { "charge", "multiplicity" } } and the invalid rows. """
        options, invalid = { }, [ ]
        for row in self.ligand_store:
            options[row[0]] = { "charge": row[4] if row[6] >= 0 else None, "multiplicity": row[5] }
            if "needs an" in row[8]: invalid.append ( "{} ({})".format ( row[0], row[8] ) )
        return options, invalid

    def _variant_model ( self, residue_name ):
        from pdynamo.opls import prep
        model = self._variant_models.get ( residue_name )
        if model is None:
            model = Gtk.ListStore ( str )
            for variant in prep.PROTONATION_CHOICES.get ( residue_name, [ ] ):
                model.append ( [ _VARIANT_TEXT.get ( variant, variant ) ] )
            self._variant_models[residue_name] = model
        return model

    def _variant_cell_data ( self, column, renderer, model, tree_iter, data = None ):
        index = model[tree_iter][_COL_INDEX]
        if self.analysis is not None and 0 <= index < len ( self.analysis.residues ):
            renderer.set_property ( "model", self._variant_model ( self.analysis.residues[index].library_name ) )

    def on_variant_edited ( self, renderer, path, text ):
        from pdynamo.opls import prep
        child_path = self.filter_model.convert_path_to_child_path ( Gtk.TreePath ( path ) )
        row = self.liststore[child_path]
        residue = self.analysis.residues[row[_COL_INDEX]]
        for variant in prep.PROTONATION_CHOICES.get ( residue.library_name, [ ] ):
            if _VARIANT_TEXT.get ( variant, variant ) == text:
                residue.variant, residue.variant_source = variant, "user"
                row[_COL_VARIANT] = text
                break

    def _row_visible ( self, model, tree_iter, data = None ):
        key = self.combo_filter.get_active_id ( ) or "all"
        if key == "all": return True
        if self.analysis is None: return True
        residue = self.analysis.residues[model[tree_iter][_COL_INDEX]]
        if key == "attention":  return model[tree_iter][_COL_ATTENTION]
        if key == "titratable":
            from pdynamo.opls import prep
            return residue.library_name in prep.PROTONATION_CHOICES
        if key == "ligands":    return residue.kind in ( "unsupported", "ligand", "other", "ion" )
        return True

    def on_combo_filter_changed ( self, widget ):
        if hasattr ( self, "filter_model" ):
            self.filter_model.refilter ( )
            self._update_count ( )

    def _update_count ( self ):
        self.label_count.set_text ( "{} of {} shown".format ( len ( self.filter_model ), len ( self.liststore ) ) )

    def _fill_table ( self ):
        from pdynamo.opls import prep
        self.liststore.clear ( )
        a = self.analysis
        breaks = { j for ( i, j, d ) in a.chain_breaks }
        for index, r in enumerate ( a.residues ):
            status, attention = [ ], False
            if r.kind == "unsupported":
                status.append ( "no OPLS template: will be parametrized (CM5)" ); attention = True
            if r.notes:
                status.append ( r.notes )
            if r.missing_heavy:
                status.append ( "missing heavy atoms: " + " ".join ( r.missing_heavy ) ); attention = True
            if r.added_atoms:
                status.append ( "placed: " + " ".join ( r.added_atoms ) )
            if r.disulfide_partner is not None:
                status.append ( "disulfide with " + a.residues[r.disulfide_partner].label )
            if index in breaks:
                status.append ( "chain break before" ); attention = True
            if r.n_hydrogens == 0 and r.kind in ( "amino acid", "water", "other" ):
                status.append ( "hydrogens will be built" )
            editable = r.library_name in prep.PROTONATION_CHOICES and r.disulfide_partner is None
            self.liststore.append ( [ index, r.label, r.kind, "; ".join ( status ) or "ok",
                                      _VARIANT_TEXT.get ( r.variant, r.variant ) if editable else "",
                                      editable, attention ] )
        self._update_count ( )
        s = a.summary ( )
        lines = [ "{} residues: {}.".format ( len ( a.residues ), ", ".join ( "{} {}".format ( n, k ) for k, n in s["kinds"].items ( ) ) ),
                  "Segments: {}.  Disulfide bridges: {}.  Chain breaks: {}.".format (
                        ", ".join ( s["segments"] ), len ( a.disulfides ), len ( a.chain_breaks ) ) ]
        if a.unsupported:
            lines.append ( "Without OPLS template (parametrized in the Ligands tab): {}.".format (
                    ", ".join ( sorted ( { a.residues[i].name for i in a.unsupported } ) ) ) )
        if s["missing_heavy"]:
            lines.append ( "Incomplete residues: {} -- tick 'Truncate' or complete them first.".format ( ", ".join ( s["missing_heavy"] ) ) )
        self.label_summary.set_text ( "\n".join ( lines ) )

    #-------------------------------------------------------------------------
    #  S I G N A L S
    #-------------------------------------------------------------------------
    def on_combobox_systems_changed ( self, widget ):
        system_id = self.combobox_systems.get_system_id ( )
        if system_id is not None:
            self.coordinates_combobox.set_model ( self.main.vobject_liststore_dict[system_id] )
            size = len ( list ( self.main.vobject_liststore_dict[system_id] ) )
            self.coordinates_combobox.set_active ( size - 1 )
            system = self.p_session.psystem[system_id]
            self.entry_system_name.set_text ( "{}_OPLS".format ( ( system.label or "system" ).replace ( " ", "_" ) ) )
            folder = getattr ( system, "e_working_folder", None ) or os.environ.get ( 'PDYNAMO3_SCRATCH' )
            if folder and os.path.isdir ( folder ):
                self.work_folder_chooser.set_folder ( folder = folder )
        self.analysis = None
        if hasattr ( self, "ligand_store" ):
            self.ligand_store.clear ( )
        if hasattr ( self, "liststore" ):
            self.liststore.clear ( )
            self.label_summary.set_text ( "Press Analyze." )

    def on_checkbox_solvate_toggled ( self, widget ):
        self.grid_solvation.set_sensitive ( widget.get_active ( ) )

    def on_button_cancel_clicked ( self, widget ):
        if not self.running: self.close_window ( )

    def _selected_system_and_coordinates ( self ):
        from pScientific.Geometry3 import Coordinates3
        system_id = self.combobox_systems.get_system_id ( )
        if system_id is None: return None, None, None
        system = self.p_session.psystem[system_id]
        vobject_id = self.coordinates_combobox.get_vobject_id ( )
        vobject = self.vm_session.vm_objects_dic.get ( vobject_id ) if vobject_id is not None else None
        coordinates3 = None
        if vobject is not None and len ( vobject.atoms ) == len ( system.atoms ):
            xyz = self.p_session.get_coordinates_from_vobject ( vobject = vobject, frame = -1 )
            coordinates3 = Coordinates3.WithExtent ( len ( xyz ) )
            for i, ( x, y, z ) in enumerate ( xyz ):
                coordinates3[i, 0], coordinates3[i, 1], coordinates3[i, 2] = x, y, z
        return system_id, system, coordinates3

    def on_button_analyze_clicked ( self, widget = None ):
        from pdynamo.opls import prep
        system_id, system, coordinates3 = self._selected_system_and_coordinates ( )
        if system is None:
            self.main.simple_dialog.info ( msg = 'Please select a system.' )
            return False
        try:
            self.analysis = prep.analyze_system ( system, coordinates3 = coordinates3 )
        except Exception as error:
            traceback.print_exc ( )
            self.analysis = None
            self.main.simple_dialog.error ( msg = 'Analysis failed:\n{}'.format ( error ) )
            return False
        self._fill_table ( )
        self._fill_ligand_table ( )
        self.label_status.set_text ( 'Analysis done. Review the residues (and ligand charges), then press Prepare.' )
        return True

    #-------------------------------------------------------------------------
    #  R U N
    #-------------------------------------------------------------------------
    def _options ( self ):
        solvation = None
        if self.checkbox_solvate.get_active ( ):
            solvation = { "shape"        : self.combo_shape.get_active_id ( ) or "cube",
                          "padding"      : self.spin_padding.get_value ( ),
                          "align"        : self.checkbox_align.get_active ( ),
                          "neutralize"   : self.checkbox_neutralize.get_active ( ),
                          "concentration": self.spin_concentration.get_value ( ),
                          "cation"       : self.combo_cation.get_active_id ( ) or "Na+",
                          "anion"        : self.combo_anion.get_active_id ( ) or "Cl-",
                          "closeness"    : self.spin_closeness.get_value ( ) }
        return { "solvation"               : solvation,
                 "parametrize_unsupported" : self.checkbox_parametrize.get_active ( ),
                 "ligand_ph"               : self.spin_ph.get_value ( ),
                 "cm5_scale"               : self.spin_cm5_scale.get_value ( ),
                 "xtb"                     : self.entry_xtb.get_text ( ).strip ( ) or None }

    def on_button_run_clicked ( self, widget ):
        from pdynamo.opls import prep
        if self.running: return
        if self.analysis is None and not self.on_button_analyze_clicked ( ):
            return
        system_id, system, _ = self._selected_system_and_coordinates ( )
        a = self.analysis
        if a.unsupported and not self.checkbox_parametrize.get_active ( ):
            self.main.simple_dialog.info ( msg = 'Residues without OPLS templates: {}.\nEnable "Parametrize" in the Ligands tab.'.format (
                    ", ".join ( sorted ( { a.residues[i].name for i in a.unsupported } ) ) ) )
            return
        incomplete = [ i for i, r in enumerate ( a.residues ) if r.missing_heavy ]
        if incomplete and not self.checkbox_truncate.get_active ( ):
            self.main.simple_dialog.info ( msg = 'Residues with missing heavy atoms: {}.\nTick "Truncate residues..." in the Residues tab or complete them first.'.format (
                    ", ".join ( a.residues[i].label for i in incomplete ) ) )
            return
        options = self._options ( )
        if a.unsupported and options["parametrize_unsupported"]:
            xtb = options["xtb"]
            if not xtb or not ( os.path.isfile ( xtb ) and os.access ( xtb, os.X_OK ) ):
                self.main.simple_dialog.info ( msg = 'Ligand charges (CM5) need the xtb executable: set its path in the Ligands tab.' )
                return
            self.vm_session.vm_config.gl_parameters['xtb_command'] = xtb
        ligand_options, invalid = self._ligand_options ( )
        if invalid:
            self.main.simple_dialog.info ( msg = 'Charge and multiplicity are incompatible for: {}.'.format ( "; ".join ( invalid ) ) )
            return
        a.ligand_options = ligand_options
        name = self.entry_system_name.get_text ( ).strip ( ) or 'system_OPLS'
        folder = self.work_folder_chooser.get_folder ( )
        if not folder or not os.path.isdir ( folder ):
            self.main.simple_dialog.info ( msg = 'Please choose a working folder (Output tab).' )
            return
        safe = "".join ( c if c.isalnum ( ) or c in "_-." else "_" for c in name )
        work_dir = os.path.join ( folder, safe + "_opls" )
        for i in incomplete:
            prep.truncate_residue ( a, i )

        self.running = True
        self.button_run.set_sensitive ( False )
        self.button_analyze.set_sensitive ( False )
        self.spinner.show ( ); self.spinner.start ( )
        self.label_status.set_text ( 'Preparing...' )

        def progress ( text ):
            GLib.idle_add ( self.label_status.set_text, text )

        def worker ( ):
            try:
                new_system, report = prep.prepare_opls_system ( system, analysis = a, work_dir = work_dir, label = name,
                                                                progress = progress, **options )
                error = None
            except Exception as exc:
                traceback.print_exc ( )
                new_system, report, error = None, None, exc
            GLib.idle_add ( self._finish, new_system, report, error, name, folder, work_dir )

        threading.Thread ( target = worker, daemon = True ).start ( )

    def _finish ( self, new_system, report, error, name, folder, work_dir ):
        from pdynamo.opls import prep
        self.running = False
        self.spinner.stop ( ); self.spinner.hide ( )
        self.button_run.set_sensitive ( True )
        self.button_analyze.set_sensitive ( True )
        analysis, self.analysis = self.analysis, None          # consumed (ligands/truncations applied)
        if error is not None:
            self.label_status.set_text ( 'Failed.' )
            self.main.simple_dialog.error_details ( parent = self.window, msg = 'The OPLS preparation failed:\n{}'.format ( error ),
                                                    details = traceback.format_exception_only ( type ( error ), error )[-1],
                                                    title = 'Prepare OPLS System' )
            return False
        text = prep.report_markdown ( report, analysis, title = "OPLS system preparation: {}".format ( name ) )
        report_path = os.path.join ( work_dir, "prepare_opls_report.md" )
        try:
            os.makedirs ( work_dir, exist_ok = True )
            with open ( report_path, "w" ) as handle: handle.write ( text )
        except OSError:
            report_path = None
        if not report["ok"]:
            self.label_status.set_text ( 'Not ready -- see the details.' )
            self.main.simple_dialog.error_details ( parent = self.window, msg = 'The system could not be fully typed/parametrized.',
                                                    details = text, title = 'Prepare OPLS System' )
            return False
        try:
            self.p_session.add_prepared_pdynamo_system ( new_system, name = name, tag = 'OPLS', working_folder = folder,
                                                         input_files = { 'opls_parameters': report.get ( "parameter_folder" ),
                                                                         'report': report_path } )
        except Exception as exc:
            traceback.print_exc ( )
            self.main.simple_dialog.error ( msg = 'The system was prepared but could not be added:\n{}'.format ( exc ) )
            return False
        solvation = report.get ( "solvation" )
        extra = ""
        if solvation:
            extra = "\n{} waters, ions: {}.".format ( solvation["waters"], ", ".join ( "{} {}".format ( n, k ) for k, n in solvation["ions"].items ( ) ) or "none" )
        if report.get ( "ligands" ):
            extra += "\nLigands parametrized: {} ({} estimated terms).".format ( ", ".join ( report["ligands"] ), len ( report.get ( "estimated_terms", [ ] ) ) )
        quality = report.get ( "quality" ) or { }
        for path, data in quality.items ( ):
            lv = data["levels"]
            extra += "\nParameter quality {}: {} green, {} yellow, {} red.".format ( path, lv["green"], lv["yellow"], lv["red"] )
            if data["torsion_scans"]:
                extra += " QC torsion scan suggested for: {}.".format ( ", ".join ( s["bond"] for s in data["torsion_scans"] ) )

        self.label_status.set_text ( 'Done: "{}" added ({} atoms).'.format ( name, report["atoms"] ) )
        if self.checkbox_open_parameters.get_active ( ) and hasattr ( self.main, "opls_parameters_window" ):
            try:
                self.main.opls_parameters_window.open_window ( e_id = getattr ( new_system, "e_id", None ) )
            except Exception:
                traceback.print_exc ( )
        self.main.simple_dialog.info ( msg = '"{}" was added: {} atoms, total charge {:+d}.{}\n\nReport: {}'.format (
                name, report["atoms"], report.get ( "total_charge", report["formal_charge"] ), extra, report_path ), title = 'Prepare OPLS System' )
        return False
