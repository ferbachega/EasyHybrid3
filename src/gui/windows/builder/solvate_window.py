#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  solvate_window.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  [EN] 2026-10-03: Builder "Solvate" window -- front end of
#  builder_solvation.solvate_builder_object() (algorithm: solvation.py).
#  Instantiated once on main_window (main.solvate_window), opened from the
#  Builder sidebar on the molecule being edited.
# ============================================================================
import os
import numpy as np
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk

from gui.windows.builder import solvation
from gui.windows.builder import solvent_library
from gui.windows.builder.builder_solvation import solvate_builder_object, solute_charge, SolvationError

_SHAPE_LABELS = ( ( "cube",   "Cube (periodic)" ),
                  ( "box",    "Rectangular box (periodic)" ),
                  ( "sphere", "Sphere (non-periodic droplet)" ) )


class SolvateWindow ( ):
    """ Builder Solvate window. """

    def __init__ ( self, main = None ):
        self.main          = main
        self.vm_session    = main.vm_session
        self.visible       = False
        self.target_object = None

    def open_window ( self, target_object ):
        self.target_object = target_object
        if not self.visible:
            self.builder = Gtk.Builder ( )
            self.builder.add_from_file ( os.path.join ( self.main.home, 'src/gui/windows/builder/solvate_window.glade' ) )
            self.builder.connect_signals ( self )
            get = self.builder.get_object
            self.window         = get ( 'window' )
            self.header_label   = get ( 'header_label' )
            self.estimate_label = get ( 'estimate_label' )
            self.status_label   = get ( 'status_label' )
            self.shape_combo    = get ( 'shape_combo' )
            self.solvent_combo  = get ( 'solvent_combo' )
            self.padding_radio  = get ( 'padding_radio' )
            self.padding_spin   = get ( 'padding_spin' )
            self.box_size_radio = get ( 'box_size_radio' )
            self.box_size_spin  = get ( 'box_size_spin' )
            self.box_b_spin     = get ( 'box_b_spin' )
            self.box_c_spin     = get ( 'box_c_spin' )
            self.align_check    = get ( 'align_check' )
            self.ion_surface_label = get ( 'ion_surface_label' )
            self.ion_surface_spin  = get ( 'ion_surface_spin' )
            self.closeness_spin = get ( 'closeness_spin' )
            self.neutralize_check   = get ( 'neutralize_check' )
            self.cation_combo       = get ( 'cation_combo' )
            self.anion_combo        = get ( 'anion_combo' )
            self.concentration_spin = get ( 'concentration_spin' )
            self.ion_solute_spin    = get ( 'ion_solute_spin' )
            self.ion_ion_spin       = get ( 'ion_ion_spin' )
            self._filling = True
            for key, label in _SHAPE_LABELS:
                self.shape_combo.append ( key, label )
            self.shape_combo.set_active_id ( "cube" )
            self._fill_solvent_combo ( )
            for key in solvation.CATIONS:
                self.cation_combo.append ( key, key )
            self.cation_combo.set_active_id ( "Na+" )
            for key in solvation.ANIONS:
                self.anion_combo.append ( key, key )
            self.anion_combo.set_active_id ( "Cl-" )
            self._filling = False
            try:
                from gui.theme import adapt_icons_for_theme
                adapt_icons_for_theme ( self.window )
            except Exception:
                pass
            self.window.show_all ( )
            self.visible = True
        else:
            self.window.present ( )
        self.status_label.set_text ( "" )
        self._refresh ( )

    def close_window ( self, *args ):
        if self.visible and getattr ( self, "window", None ) is not None:
            self.window.destroy ( )
        self.visible = False
        return True

    def on_close_button_clicked ( self, *args ):
        return self.close_window ( )

    # ------------------------------------------------------------------
    def _fill_solvent_combo ( self, select_path = None ):
        """ Library solvents (solvent_library.discover_solvents()); id = path. """
        previous = select_path or self.solvent_combo.get_active_id ( )
        self.solvent_combo.remove_all ( )
        self._solvent_paths = [ ]
        for key, label, path, kind in solvent_library.discover_solvents ( ):
            suffix = "" if kind == "box" else "  (box built from the molecule)"
            self.solvent_combo.append ( path, label + suffix )
            self._solvent_paths.append ( path )
        if previous and previous in self._solvent_paths:
            self.solvent_combo.set_active_id ( previous )
        elif self._solvent_paths:
            self.solvent_combo.set_active ( 0 )

    def _solvent_info ( self ):
        path = self.solvent_combo.get_active_id ( )
        if not path:
            return None, None
        try:
            return path, solvent_library.quick_info ( path )
        except solvent_library.SolventError as error:
            return path, { "error": str ( error ) }

    def on_solvent_changed ( self, combo ):
        if not getattr ( self, "_filling", False ) and self.visible:
            self._refresh ( )

    def _parameters ( self ):
        shape = self.shape_combo.get_active_id ( ) or "cube"
        box_size = None
        if self.box_size_radio.get_active ( ):
            if shape == "box":
                box_size = ( self.box_size_spin.get_value ( ), self.box_b_spin.get_value ( ),
                             self.box_c_spin.get_value ( ) )
            else:
                box_size = self.box_size_spin.get_value ( )
        return dict (
            shape               = shape,
            padding             = self.padding_spin.get_value ( ),
            box_size            = box_size,
            align               = shape == "box" and self.align_check.get_active ( ),
            ion_surface_distance = self.ion_surface_spin.get_value ( ),
            closeness           = self.closeness_spin.get_value ( ),
            neutralize          = self.neutralize_check.get_active ( ),
            concentration       = self.concentration_spin.get_value ( ) / 1000.0,
            cation              = self.cation_combo.get_active_id ( ) or "Na+",
            anion               = self.anion_combo.get_active_id ( ) or "Cl-",
            ion_solute_distance = self.ion_solute_spin.get_value ( ),
            ion_ion_distance    = self.ion_ion_spin.get_value ( ),
        )

    def _refresh ( self, warn = True ):
        """ Header + a quick estimate (box edge, approximate number of
            waters and ions) without running the real placement. """
        obj = self.target_object
        if obj is None or not getattr ( obj, "atoms", None ):
            self.header_label.set_markup ( "<b>No molecule in the Builder.</b>" )
            self.estimate_label.set_text ( "" )
            return
        charge = solute_charge ( obj )
        self.header_label.set_markup ( "<b>{}</b>: {} atoms, total formal charge {:+d}".format (
                obj.name, len ( obj.atoms ), charge ) )
        solvation_info = getattr ( obj, "builder_solvation", None )
        cell = getattr ( obj, "builder_cell", None )
        if cell or solvation_info:
            if cell:
                text = "Current periodic cell: {:.2f} x {:.2f} x {:.2f} \u00c5.".format ( *cell[:3] )
            else:
                text = "Solvated in a non-periodic sphere of radius {:.2f} \u00c5.".format (
                        solvation_info.get ( "radius", 0.0 ) )
            self.estimate_label.set_text ( text )
            if warn:
                self.status_label.set_text ( "This molecule was already solvated; undo the previous "
                                             "solvation before solvating again." )
            self.builder.get_object ( "solvate_button" ).set_sensitive ( False )
            return
        self.builder.get_object ( "solvate_button" ).set_sensitive ( True )

        shape = self.shape_combo.get_active_id ( ) or "cube"
        fixed = self.box_size_radio.get_active ( )
        self.padding_spin.set_sensitive ( not fixed )
        self.box_size_spin.set_sensitive ( fixed )
        self.box_b_spin.set_visible ( shape == "box" )
        self.box_c_spin.set_visible ( shape == "box" )
        self.box_b_spin.set_sensitive ( fixed )
        self.box_c_spin.set_sensitive ( fixed )
        self.box_size_radio.set_label ( { "cube": "Box edge (\u00c5):", "box": "Edges a, b, c (\u00c5):",
                                          "sphere": "Radius (\u00c5):" }[shape] )
        self.padding_radio.set_label ( "Margin around solute (\u00c5):" )
        self.align_check.set_visible ( shape == "box" )
        self.ion_surface_label.set_visible ( shape == "sphere" )
        self.ion_surface_spin.set_visible ( shape == "sphere" )
        if warn:
            self.status_label.set_text ( "" )
        try:
            p = self._parameters ( )
            xyz = obj.frames[0][sorted ( obj.atoms.keys ( ) )]
            region = solvation.region_for_solute ( xyz, p["shape"], p["padding"], p["box_size"], p["align"] )
            # 0.0334 waters/A3 at 1 g/cm3; periodic boxes lose ~6 % at the
            # seams; each solute heavy atom (+ the 2.4 A exclusion shell)
            # removes ~25 A3 -- calibrated on ubiquitin (box, sphere, cube
            # all within ~6 % of the real count)
            path, info = self._solvent_info ( )
            if info is None or "error" in info or not info.get ( "number_density" ):
                self.estimate_label.set_text ( ( info or { } ).get ( "error", "Choose a solvent." ) )
                return
            # number density of the solvent; periodic boxes lose ~6 % at the
            # seams; each solute heavy atom (+ the exclusion shell) removes a
            # volume growing with the solvent molecule's size, ~k x V_mol^(1/3)
            # (k = 8 periodic, 6.2 sphere) -- fitted on ubiquitin with the 13
            # library solvents, estimates within ~10 % of the real count
            nd      = info["number_density"]
            density = nd * ( 0.94 if region.periodic else 0.99 )
            excluded_per_heavy = ( 8.0 if region.periodic else 6.2 ) * ( 1.0 / nd ) ** ( 1.0 / 3.0 )
            heavy = sum ( 1 for a in obj.atoms.values ( ) if a.symbol != "H" )
            n_mol = max ( 0, int ( density * ( region.volume - excluded_per_heavy * heavy ) ) )
            molarity = info["density"] * 1000.0 / info["molar_mass"]
            n_cat, n_an = solvation.ion_counts ( charge, n_mol, p["concentration"], p["neutralize"], molarity )
            what = "waters" if info["n_atoms"] == 3 and info["residue"] in ( "HOH", "WAT", "SOL", "TIP3" ) \
                   else "{} molecules".format ( info["residue"] )
            text = "Estimate: {}, ~{} {}, {} {} + {} {} (~{} atoms in total){}.".format (
                    region.describe ( ), n_mol, what, n_cat, p["cation"], n_an, p["anion"],
                    len ( obj.atoms ) + info["n_atoms"] * n_mol,
                    "" if region.periodic else "; no periodic cell" )
            if not info["cached"]:
                text += " The {} box is built the first time it is used (a few seconds).".format ( info["name"] )
            self.estimate_label.set_text ( text )
        except SolvationError as error:
            self.estimate_label.set_text ( str ( error ) )

    def on_parameter_changed ( self, *args ):
        if not getattr ( self, "_filling", False ) and self.visible:
            self._refresh ( )

    def on_solvate_button_clicked ( self, button ):
        obj = self.target_object
        if obj is None or not getattr ( obj, "atoms", None ):
            self.status_label.set_text ( "There is no molecule in the Builder to solvate." )
            return
        path = self.solvent_combo.get_active_id ( )
        if not path:
            self.status_label.set_text ( "Choose a solvent." )
            return
        def pump ( ):
            while Gtk.events_pending ( ):
                Gtk.main_iteration ( )
        def progress ( fraction ):
            self.status_label.set_text ( "Building the solvent box (first use only, cached): {:.0f} %".format (
                    100.0 * min ( 1.0, max ( 0.0, fraction ) ) ) )
            pump ( )
        self.status_label.set_text ( "Loading the solvent..." )
        pump ( )
        try:
            solvent = solvent_library.load_solvent ( path, progress = progress )
        except solvent_library.SolventError as error:
            self.status_label.set_text ( str ( error ) )
            return
        self.status_label.set_text ( "Solvating..." )
        pump ( )
        try:
            result = solvate_builder_object ( obj, solvent = solvent, **self._parameters ( ) )
        except SolvationError as error:
            self.status_label.set_text ( str ( error ) )
            return
        except Exception as error:
            import traceback
            traceback.print_exc ( )
            self.status_label.set_text ( "Solvation failed: {}".format ( error ) )
            return
        message = "Solvated: " + result.summary ( ) + ". Minimise and equilibrate (NPT) before production."
        self.status_label.set_text ( message )
        # look at the whole box / sphere
        try:
            if result.cell is not None:
                centre = np.array ( result.cell[:3], dtype = np.float32 ) / 2.0
            else:
                centre = np.array ( result.sphere[0], dtype = np.float32 )
            self.vm_session.vm_glcore.center_on_coordinates ( obj, centre )
        except Exception:
            pass
        statusbar = getattr ( self.main, "statusbar_main", None )
        if statusbar is not None:
            statusbar.push ( 1, message )
        self._refresh ( warn = False )

    # ------------------------------------------------------------------
    def on_add_solvent_button_clicked ( self, button ):
        """ Adds a .mol2 file to the solvent library (solvent_library.
            add_solvent_file()): a periodic box with a CRYSIN section is
            copied as is; for a single molecule the name, residue name and
            density are asked for. """
        chooser = Gtk.FileChooserDialog ( title = "Add a solvent (.mol2)", parent = self.window,
                                          action = Gtk.FileChooserAction.OPEN )
        chooser.add_buttons ( Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK )
        mol2_filter = Gtk.FileFilter ( )
        mol2_filter.set_name ( "Tripos MOL2 (*.mol2)" )
        mol2_filter.add_pattern ( "*.mol2" )
        mol2_filter.add_pattern ( "*.MOL2" )
        chooser.add_filter ( mol2_filter )
        response = chooser.run ( )
        source = chooser.get_filename ( )
        chooser.destroy ( )
        if response != Gtk.ResponseType.OK or not source:
            return
        try:
            info = solvent_library.quick_info ( source )
        except Exception as error:
            self.status_label.set_text ( "Could not read {}: {}".format ( os.path.basename ( source ), error ) )
            return
        values = self._ask_solvent_properties ( info )
        if values is None:
            return
        name, residue, density = values
        if info["kind"] == "molecule" and not density:
            self.status_label.set_text ( "A single-molecule solvent needs its density." )
            return
        try:
            target = solvent_library.add_solvent_file ( source, name = name, residue = residue, density = density )
            solvent_library.quick_info ( target )
        except Exception as error:
            self.status_label.set_text ( "Could not add the solvent: {}".format ( error ) )
            return
        self._filling = True
        self._fill_solvent_combo ( select_path = target )
        self._filling = False
        self.status_label.set_text ( "Added \"{}\" to the solvent library ({}).".format (
                name or info["name"], os.path.basename ( target ) ) )
        self._refresh ( warn = False )

    def _ask_solvent_properties ( self, info ):
        """ Small dialog: name, residue name, density (g/cm3). Returns
            (name, residue, density or None) or None if cancelled. """
        dialog = Gtk.Dialog ( title = "New solvent", transient_for = self.window, modal = True )
        dialog.add_buttons ( Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OK, Gtk.ResponseType.OK )
        grid = Gtk.Grid ( row_spacing = 6, column_spacing = 8, margin = 10 )
        kind = "pre-equilibrated periodic box" if info["kind"] == "box" else "single molecule (a box will be built at this density)"
        grid.attach ( Gtk.Label ( label = "File type: " + kind, xalign = 0 ), 0, 0, 2, 1 )
        entry_name = Gtk.Entry ( text = info["name"] or "" )
        entry_residue = Gtk.Entry ( text = info["residue"] or "SOL", max_length = 4 )
        density_spin = Gtk.SpinButton.new_with_range ( 0.05, 5.0, 0.001 )
        density_spin.set_digits ( 3 )
        density_spin.set_value ( info["density"] or 1.0 )
        for row, ( label, widget ) in enumerate ( ( ( "Name:", entry_name ), ( "Residue name:", entry_residue ),
                                                    ( "Density (g/cm\u00b3):", density_spin ) ), start = 1 ):
            grid.attach ( Gtk.Label ( label = label, xalign = 0 ), 0, row, 1, 1 )
            grid.attach ( widget, 1, row, 1, 1 )
        if info["kind"] == "box" and info["density"]:
            density_spin.set_sensitive ( False )      # the box defines its own density
        dialog.get_content_area ( ).add ( grid )
        dialog.show_all ( )
        response = dialog.run ( )
        values = ( entry_name.get_text ( ).strip ( ) or None,
                   entry_residue.get_text ( ).strip ( ).upper ( ) or None,
                   None if info["kind"] == "box" and info["density"] else density_spin.get_value ( ) )
        dialog.destroy ( )
        return values if response == Gtk.ResponseType.OK else None
