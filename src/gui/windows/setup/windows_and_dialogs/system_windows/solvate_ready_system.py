#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: "Solvate System" window
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description (2026-10-21, user's request): System > Solvate System...
#  solvates a system that is already prepared (OPLS or DYFF force field) and
#  adds the result as a new system. Replaces the old Edit > Solvate System /
#  Make Solvent Box windows (system_windows/solvate_system.py, which needed a
#  pre-built solvent-box system and ion .pkl files). The work is done by
#  pdynamo/solvate_system.py in a worker thread; the preview (box, waters,
#  ions) uses only the fast placement step.
#

import copy
import os
import threading
import traceback

import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk, GLib

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox


class SolvateReadySystemWindow:

    def __init__ ( self, main = None ):
        self.main       = main
        self.home       = main.home
        self.p_session  = main.p_session
        self.vm_session = main.vm_session
        self.Visible    = False
        self._busy      = False
        self._preview_id = None

    # ------------------------------------------------------------------
    def open_window ( self ):
        if self.Visible:
            self.window.present ( ); return
        self.builder = Gtk.Builder ( )
        self.builder.add_from_file ( os.path.join ( self.home, "src/gui/windows/setup/solvate_ready_system_window.glade" ) )
        self.builder.connect_signals ( self )
        get = self.builder.get_object
        self.window = get ( "window" )
        self.window.set_transient_for ( self.main.window )
        self.window.connect ( "delete-event", lambda w, e: ( self.close_window ( ), True )[1] )
        self.combo_shape, self.combo_cation, self.combo_anion = get ( "combo_shape" ), get ( "combo_cation" ), get ( "combo_anion" )
        self.spin_padding, self.spin_closeness = get ( "spin_padding" ), get ( "spin_closeness" )
        self.spin_concentration, self.check_neutralize, self.check_align = get ( "spin_concentration" ), get ( "check_neutralize" ), get ( "check_align" )
        self.entry_name, self.label_preview = get ( "entry_name" ), get ( "label_preview" )
        self.label_status, self.label_system_info = get ( "label_status" ), get ( "label_system_info" )
        self.button_solvate = get ( "button_solvate" )

        self.combobox_systems = SystemComboBox ( self.main )
        get ( "box_system" ).pack_start ( self.combobox_systems, True, True, 0 )
        self.coordinates_combobox = CoordinatesComboBox ( )
        get ( "box_coordinates" ).pack_start ( self.coordinates_combobox, True, True, 0 )
        self.combobox_systems.connect ( "changed", self.on_system_changed )
        self.coordinates_combobox.connect ( "changed", self.on_options_changed )

        self.window.show_all ( )
        self.Visible = True
        if self.p_session.psystem:
            self.combobox_systems.set_active_system ( e_id = self.p_session.active_id )
        self.on_system_changed ( )

    def close_window ( self, *args ):
        if not self.Visible: return
        if self._busy:
            self.label_status.set_text ( "Solvating... the window closes when it finishes." ); return
        self.window.destroy ( )
        self.Visible = False

    # ------------------------------------------------------------------
    def _system ( self ):
        e_id = self.combobox_systems.get_system_id ( )
        return self.p_session.psystem.get ( e_id ) if e_id is not None else None

    def on_system_changed ( self, *args ):
        from pdynamo import solvate_system as SS
        system = self._system ( )
        if system is None:
            self.label_system_info.set_text ( "" ); return
        e_id = system.e_id
        if e_id in self.main.vobject_liststore_dict:
            model = self.main.vobject_liststore_dict[e_id]
            self.coordinates_combobox.set_model ( model )
            self.coordinates_combobox.set_active ( len ( model ) - 1 )
        self.entry_name.set_text ( "{}_solvated".format ( system.label ) )
        kind = SS.force_field_kind ( system )
        for combo, ions in ( ( self.combo_cation, 0 ), ( self.combo_anion, 1 ) ):
            combo.remove_all ( )
            for ion in SS.IONS_BY_FORCE_FIELD.get ( kind or "OPLS" )[ions]:
                combo.append ( ion, ion )
            combo.set_active ( 0 )
        try:
            SS.check ( system )
            charge, exact = SS.solute_charge ( system )
            text = "{} atoms · {} · charge {:+.3f}".format ( len ( system.atoms ), kind, exact )
            if abs ( exact - charge ) > 0.01:
                text += "  (not an integer: the counter-ions neutralize {:+d})".format ( charge )
            self.label_system_info.set_text ( text )
            self.button_solvate.set_sensitive ( True )
        except SS.SolvateSystemError as error:
            self.label_system_info.set_markup ( '<span foreground="#e5534b">{}</span>'.format ( GLib.markup_escape_text ( str ( error ) ) ) )
            self.button_solvate.set_sensitive ( False )
        self.on_options_changed ( )

    def _options ( self ):
        shape = self.combo_shape.get_active_id ( ) or "cube"
        return { "shape": shape, "padding": self.spin_padding.get_value ( ), "closeness": self.spin_closeness.get_value ( ),
                 "neutralize": self.check_neutralize.get_active ( ), "concentration": self.spin_concentration.get_value ( ),
                 "cation": self.combo_cation.get_active_id ( ) or "Na+", "anion": self.combo_anion.get_active_id ( ) or "Cl-",
                 "align": self.check_align.get_active ( ) and shape == "box", "seed": 1 }

    def _coordinates ( self, system ):
        """ Coordinates3 of the chosen object/frame (the system itself is not changed). """
        from pScientific.Geometry3 import Coordinates3
        vobject_id = self.coordinates_combobox.get_vobject_id ( )
        vobject = self.vm_session.vm_objects_dic.get ( vobject_id ) if vobject_id is not None else None
        if vobject is None or len ( vobject.atoms ) != len ( system.atoms ):
            return system.coordinates3
        frame = min ( getattr ( self.vm_session, "frame", 0 ) or 0, len ( vobject.frames ) - 1 )
        xyz = vobject.frames[frame]
        crd = Coordinates3.WithExtent ( len ( system.atoms ) )
        for i in range ( len ( system.atoms ) ):
            crd[i, 0], crd[i, 1], crd[i, 2] = float ( xyz[i][0] ), float ( xyz[i][1] ), float ( xyz[i][2] )
        return crd

    def on_options_changed ( self, *args ):
        self.check_align.set_sensitive ( ( self.combo_shape.get_active_id ( ) or "cube" ) == "box" )
        # . the preview runs after the user stops changing values
        if self._preview_id is not None:
            GLib.source_remove ( self._preview_id )
        self._preview_id = GLib.timeout_add ( 300, self._update_preview )

    def _update_preview ( self ):
        from pdynamo import solvate_system as SS
        self._preview_id = None
        if not self.Visible: return False
        system = self._system ( )
        if system is None or getattr ( system, "mmModel", None ) is None or SS.force_field_kind ( system ) is None:
            self.label_preview.set_text ( "" ); return False
        try:
            options = self._options ( )
            charge, exact = SS.solute_charge ( system )
            result = SS.place ( SS.system_xyz ( system, self._coordinates ( system ) ), charge, **options )
        except Exception as error:
            self.label_preview.set_text ( "Preview not available: {}".format ( error ) ); return False
        ions = { }
        for ( ion, *_ ) in result.ions: ions[ion] = ions.get ( ion, 0 ) + 1
        if result.cell:
            box = "Box {:.1f} x {:.1f} x {:.1f} Å".format ( *result.cell[:3] )
        else:
            box = "Droplet radius {:.1f} Å".format ( result.sphere[1] )
        atoms = len ( system.atoms ) + 3 * result.n_waters + len ( result.ions )
        self.label_preview.set_markup ( "<b>Preview:</b> {} · {} waters · ions: {} · <b>{} atoms</b>".format (
            box, result.n_waters, ", ".join ( "{} {}".format ( n, k ) for k, n in sorted ( ions.items ( ) ) ) or "none", atoms ) )
        return False

    # ------------------------------------------------------------------
    def on_button_solvate_clicked ( self, button ):
        from pdynamo import solvate_system as SS
        if self._busy: return
        system = self._system ( )
        try:
            SS.check ( system )
        except SS.SolvateSystemError as error:
            self.main.simple_dialog.error ( msg = str ( error ) ); return
        name = self.entry_name.get_text ( ).strip ( ) or "{}_solvated".format ( system.label )
        options = self._options ( )
        coordinates = self._coordinates ( system )
        self._busy = True
        self.button_solvate.set_sensitive ( False )
        self.label_status.set_text ( "Solvating..." )
        self.main.statusbar_main.push ( 1, "Solvate System: {}...".format ( system.label ) )

        def say ( text ):
            GLib.idle_add ( lambda: ( self.Visible and self.label_status.set_text ( "Solvating: {}...".format ( text ) ), False )[1] )

        def work ( ):
            try:
                new, report = SS.solvate ( system, coordinates3 = coordinates, label = name, progress = say, **options )
                GLib.idle_add ( self._done, system, new, report, name, None )
            except SS.SolvateSystemError as error:
                GLib.idle_add ( self._done, system, None, None, name, str ( error ) )
            except Exception as error:
                traceback.print_exc ( )
                GLib.idle_add ( self._done, system, None, None, name, "{}: {}".format ( type ( error ).__name__, error ) )
        threading.Thread ( target = work, daemon = True ).start ( )

    def _done ( self, source, new, report, name, error ):
        from pdynamo import solvate_system as SS
        self._busy = False
        if self.Visible:
            self.button_solvate.set_sensitive ( True )
        if error:
            if self.Visible: self.label_status.set_text ( "Failed." )
            self.main.simple_dialog.error ( msg = "Solvate System failed:\n{}".format ( error ) )
            return False
        input_files = dict ( getattr ( source, "e_input_files", None ) or { } )
        system = self.p_session.add_prepared_pdynamo_system ( new, name = name, tag = "solvated",
                                                              working_folder = getattr ( source, "e_working_folder", None ),
                                                              input_files = input_files )
        notes = self._carry_over ( source, system )
        text = SS.report_text ( report )
        if notes: text += "\n" + "\n".join ( notes )
        text += "\nMinimize the new system (Simulate > Geometry Optimization) before running MD."
        self.main.bottom_notebook.status_teeview_add_new_item (
            message = "Solvated '{}': {} waters, ions {}, {} atoms.".format ( source.label, report["waters"],
                      ", ".join ( "{} {}".format ( n, k ) for k, n in sorted ( report["ions"].items ( ) ) ) or "none", report["atoms"] ),
            system = system )
        self.main.refresh_main_statusbar ( )
        if self.Visible:
            self.label_status.set_text ( "Done: new system '{}'.".format ( system.label ) )
        self.main.simple_dialog.info ( msg = text, title = "Solvate System" )
        return False

    def _carry_over ( self, source, system ):
        """ Session data of the solute that stays valid (the solute atoms keep their
            indices): selections, fixed atoms and, for DYFF, the data kept on the
            object that lets the MM model be rebuilt the same way. """
        notes = [ ]
        if getattr ( source, "e_selections", None ):
            system.e_selections = copy.deepcopy ( source.e_selections )
            notes.append ( "Selections copied: {}.".format ( ", ".join ( source.e_selections ) ) )
        vobject = next ( ( v for v in self.vm_session.vm_objects_dic.values ( ) if getattr ( v, "e_id", None ) == system.e_id ), None )
        fixed = list ( getattr ( source, "e_fixed_table", None ) or [ ] )
        if fixed:
            from pCore import Selection
            selection = Selection.FromIterable ( fixed )
            system.e_fixed_table = list ( selection )
            system.freeAtoms = selection.Complement ( upperBound = len ( system.atoms ) )
            if vobject is not None:
                try:
                    self.p_session._apply_fixed_representation_to_vobject ( vismol_object = vobject )
                except Exception:
                    traceback.print_exc ( )
            notes.append ( "{} fixed atoms kept fixed (the solvent is free).".format ( len ( fixed ) ) )
        if getattr ( source, "e_restraints_dict", None ):
            notes.append ( "Restraints were not copied: define them again on the new system." )
        if vobject is not None and str ( getattr ( system.mmModel, "forceField", "" ) ).upper ( ) == "DYFF":
            self._carry_dyff_object_data ( source, system, vobject )
        if vobject is not None:
            self.vm_session.vm_glcore.queue_draw ( )
        return notes

    def _carry_dyff_object_data ( self, source, system, vobject ):
        from pdynamo.dyff_edits import term_signature
        old = next ( ( v for v in self.vm_session.vm_objects_dic.values ( )
                       if getattr ( v, "e_id", None ) == source.e_id and getattr ( v, "dyff_parameter_edits", None ) is not None ), None )
        old = old or next ( ( v for v in self.vm_session.vm_objects_dic.values ( ) if getattr ( v, "e_id", None ) == source.e_id ), None )
        if old is not None:
            edits = getattr ( old, "dyff_parameter_edits", None )
            if edits:
                edits = copy.deepcopy ( edits )
                if "term_signature" in edits:
                    edits["term_signature"] = term_signature ( len ( system.atoms ), list ( system.connectivity.bondIndices ) )
                vobject.dyff_parameter_edits = edits
            for attribute in ( "dyff_type_overrides", "dyff_library_options" ):
                value = getattr ( old, attribute, None )
                if value is not None:
                    setattr ( vobject, attribute, copy.deepcopy ( value ) )
        # . the charges the next DYFF definition reads from the object (solute + TIP3P + ions)
        charges = system.mmState.charges
        for i, atom in vobject.atoms.items ( ):
            if i < len ( system.atoms ):
                atom.mm_charge = float ( charges[i] )
