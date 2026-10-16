#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: DYFF Parameters window of the Builder (replaces the Atom Types window)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Tabs: Atoms (DYFF type -- editable, validated --, perceived type,
#      hybridization, formal charge, partial charge, quality), Bonds,
#      Angles, Dihedrals, Impropers and Lennard-Jones (the parameters DYFF
#      builds from the types, with the strain of each term at the current
#      geometry as the quality colour; see dyff_parameters.py). Everything
#      the Atom Types window did is here: type overrides, formal charges,
#      Protonate/Deprotonate, Clear Override, Rename Residue, the numbered
#      atom labels in the 3D view. New: partial charges (by hand or CM5 with
#      xTB) and the projection of any row on the structure with the picking
#      markers pk1-pk4 (previous/next term of the same parameter).
#      Dihedrals: "Fit with xTB scan..." scans the dihedral shown with xTB
#      and fits the row's parameter to it (torsion_fit_window.py).
#

import os
import traceback
import types as _types

import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk, GLib

from gui.windows.builder import dyff_parameters as D


_TERM_TABS = [
    ( "bond",       "Bonds",         [ ( "Types and bond tag", "types" ), ( "b0 (A)", "v0" ), ( "k (kJ/mol/A2)", "v1" ),
                                        ( "Terms", "count" ), ( "Example", "example" ), ( "Measured", "measured" ) ] ),
    ( "angle",      "Angles",        [ ( "Types and bond tags", "types" ), ( "minimum (deg)", "v0" ), ( "c n (kJ/mol)", "v1" ),
                                        ( "Terms", "count" ), ( "Example", "example" ), ( "Measured", "measured" ) ] ),
    ( "dihedral",   "Dihedrals",     [ ( "Types", "types" ), ( "c n (kJ/mol)", "v0" ),
                                        ( "Terms", "count" ), ( "Example", "example" ), ( "Measured", "measured" ) ] ),
    ( "outofplane", "Impropers",     [ ( "Types", "types" ), ( "c n (kJ/mol)", "v0" ),
                                        ( "Terms", "count" ), ( "Example", "example" ), ( "Measured", "measured" ) ] ),
    ( "lj",         "Lennard-Jones", [ ( "Type", "types" ), ( "epsilon", "v0" ), ( "sigma", "v1" ), ( "Atoms", "count" ) ] ),
]
# . Atoms tab: ( title, field, editable )
_ATOM_COLUMNS = [ ( "Atom", "atom_id", False ), ( "Symbol", "symbol", False ), ( "Name", "name", False ),
                  ( "Residue", "residue", False ), ( "DYFF type", "type", True ), ( "Perceived", "perceived", False ),
                  ( "Hybridization", "hybridization", False ), ( "Formal charge", "formal", True ),
                  ( "Partial charge", "charge", True ) ]
_FILTERS = [ ( "selected", "Selected atoms" ), ( "attention", "Needs attention (yellow/red)" ), ( "all", "All" ) ]


class DYFFParametersWindow:
    """ DYFF Parameters window (Builder). """

    def __init__ ( self, main = None ):
        self.main       = main
        self.vm_session = main.vm_session
        self.visible    = False
        self.Visible    = False
        self.target_object = None
        self.selected_atoms = [ ]
        self.tables     = None
        self.current    = None

    #-------------------------------------------------------------------------
    def open_window ( self, target_object, atom_ids = None ):
        self.target_object = target_object
        self.selected_atoms = [ target_object.atoms[i] for i in sorted ( atom_ids or [ ] ) if i in target_object.atoms ]
        if not self.visible:
            self.builder = Gtk.Builder ( )
            self.builder.add_from_file ( os.path.join ( self.main.home, 'src/gui/windows/builder/dyff_parameters_window.glade' ) )
            self.builder.connect_signals ( self )
            get = self.builder.get_object
            self.window         = get ( 'window' )
            self.notebook       = get ( 'notebook' )
            self.combo_filter   = get ( 'combo_filter' )
            self.label_charge   = get ( 'label_charge' )
            self.label_status   = get ( 'label_status' )
            self.label_legend   = get ( 'label_legend' )
            self.label_instance = get ( 'label_instance' )
            self.header_label   = get ( 'header_label' )
            self.checkbox_project = get ( 'checkbox_project' )
            self.checkbox_center  = get ( 'checkbox_center' )
            self.button_fit_torsion = get ( 'button_fit_torsion' )
            self.label_legend.set_markup ( "  ".join ( '<span background="{}">  </span> {}'.format ( D.LEVEL_COLORS[l], t )
                                                      for l, t in ( ( "green", "ok" ), ( "yellow", "check" ),
                                                                    ( "red", "strained / untyped" ), ( "blue", "manual type" ) ) ) )
            for key, text in _FILTERS: self.combo_filter.append ( key, text )
            self._build_tabs ( )
            self._build_atom_actions ( get ( 'box_atom_actions' ) )
            self.window.connect ( 'destroy', self.close_window )
            self.window.show_all ( )
            self.builder.get_object ( 'button_apply' ).hide ( )      # DYFF edits apply immediately
            self.visible = self.Visible = True
        else:
            self.window.present ( )
        self.combo_filter.set_active_id ( "selected" if self.selected_atoms else "all" )
        self.refresh ( )

    def close_window ( self, *args ):
        if not self.visible: return
        if getattr ( self, "torsion_fit_window", None ) is not None:
            self.torsion_fit_window.close_window ( )
        self.vm_session.builder_atom_types_labeled_atoms = [ ]
        self.window.destroy ( )
        self.visible = self.Visible = False
        if getattr ( self.vm_session, "vm_glcore", None ) is not None:
            self.vm_session.vm_glcore.queue_draw ( )

    def on_button_close_clicked ( self, widget ):
        self.close_window ( )

    def on_button_revert_clicked ( self, widget ):
        self.refresh ( )
        self.label_status.set_text ( "Refreshed." )

    def on_button_apply_clicked ( self, widget ):
        pass                                   # DYFF edits are applied immediately

    #-------------------------------------------------------------------------
    #  T A B S
    #-------------------------------------------------------------------------
    def _new_tree ( self, term, title, columns, editable_fields = ( ) ):
        n = len ( columns )
        store = Gtk.ListStore ( *( [ int ] + [ str ] * n + [ str, str, str, bool ] ) )
        model = store.filter_new ( )
        model.set_visible_func ( self._row_visible )
        tree = Gtk.TreeView ( model = model )
        if term == "atoms": tree.get_selection ( ).set_mode ( Gtk.SelectionMode.MULTIPLE )
        for k, ( name, field ) in enumerate ( columns ):
            renderer = Gtk.CellRendererText ( )
            if field in editable_fields:
                renderer.set_property ( "editable", True )
                renderer.set_property ( "weight", 600 )
                if term == "atoms": renderer.connect ( "edited", self.on_atom_cell_edited, field, k + 1 )
                else:               renderer.connect ( "edited", self.on_term_cell_edited, term, field )
            column = Gtk.TreeViewColumn ( name, renderer, text = k + 1 )
            column.set_resizable ( True )
            column.set_sort_column_id ( k + 1 )
            tree.append_column ( column )
        renderer = Gtk.CellRendererText ( )
        renderer.set_property ( "foreground", "#1e1e1e" )
        column = Gtk.TreeViewColumn ( "Quality", renderer, text = n + 1, background = n + 2 )
        column.set_sort_column_id ( n + 3 )
        tree.append_column ( column )
        tree.get_selection ( ).connect ( "changed", self.on_row_selected, term )
        tree.connect ( "row-activated", self.on_row_activated )
        scrolled = Gtk.ScrolledWindow ( )
        scrolled.add ( tree )
        label = Gtk.Label ( label = title )
        self.notebook.append_page ( scrolled, label )
        self.stores[term], self.filters[term], self.trees[term] = store, model, tree
        self.columns[term] = columns
        self.tab_labels[term] = ( label, title )

    def _build_tabs ( self ):
        self.stores, self.filters, self.trees, self.columns, self.tab_labels = { }, { }, { }, { }, { }
        self._new_tree ( "atoms", "Atoms", [ ( t, f ) for t, f, e in _ATOM_COLUMNS ],
                         editable_fields = [ f for t, f, e in _ATOM_COLUMNS if e ] )
        editable = { "bond": ( "v0", "v1" ), "angle": ( "v1", ), "dihedral": ( "v0", ), "outofplane": ( "v0", ) }
        for term, title, columns in _TERM_TABS:
            self._new_tree ( term, title, columns, editable_fields = editable.get ( term, ( ) ) )
        self.notebook.connect ( "switch-page", lambda nb, page, n: self._update_labeled_atoms ( n ) )

    def _build_atom_actions ( self, box ):
        for label, handler, tip in (
                ( "Protonate", self.on_protonate, "Add a proton to the selected rows' atoms (+1)" ),
                ( "Deprotonate", self.on_deprotonate, "Remove a proton from the selected rows' atoms (-1)" ),
                ( "Clear type", self.on_clear_override, "Back to the perceived DYFF type for the selected rows" ),
                ( "Rename Residue...", self.on_rename_residue, "Rename the residue(s) of the selected rows" ),
                ( "CM5 charges...", self.on_cm5_charges, "Partial charges of the whole molecule: GFN1-xTB CM5" ),
                ( "Clear charges", self.on_clear_charges, "Remove every partial charge (DYFF: all zero)" ) ):
            button = Gtk.Button ( label = label )
            button.set_tooltip_text ( tip )
            button.connect ( "clicked", handler )
            box.pack_start ( button, False, False, 0 )
        box.show_all ( )

    def _row_visible ( self, model, tree_iter, data = None ):
        key = self.combo_filter.get_active_id ( ) or "all"
        n = model.get_n_columns ( )
        level, involved = model[tree_iter][n - 2], model[tree_iter][n - 1]
        if key == "selected":  return involved
        if key == "attention": return level in ( "yellow", "red" )
        return True

    def on_combo_filter_changed ( self, widget ):
        if hasattr ( self, "filters" ):
            for model in self.filters.values ( ): model.refilter ( )
            self._update_tab_titles ( )
            self._update_labeled_atoms ( )

    def _update_tab_titles ( self ):
        for term, ( label, title ) in self.tab_labels.items ( ):
            store = self.stores[term]
            n = store.get_n_columns ( )
            red    = sum ( 1 for row in store if row[n - 2] == "red" )
            yellow = sum ( 1 for row in store if row[n - 2] == "yellow" )
            shown, total = len ( self.filters[term] ), len ( store )
            markup = "{} ({})".format ( title, total ) if shown == total else "{} ({}/{})".format ( title, shown, total )
            if term not in ( "atoms", "lj" ) and self.tables is not None:
                n_terms = sum ( r["count"] for r in self.tables.get ( term, [ ] ) )
                if n_terms != total: markup += " <small>{} terms</small>".format ( n_terms )
            if red:    markup += ' <span foreground="#e5534b">●{}</span>'.format ( red )
            if yellow: markup += ' <span foreground="#c9a227">●{}</span>'.format ( yellow )
            label.set_markup ( markup )

    #-------------------------------------------------------------------------
    #  F I L L
    #-------------------------------------------------------------------------
    def refresh ( self ):
        if self.target_object is None: return
        self.current = None
        self.button_fit_torsion.set_sensitive ( False )
        self.selected_atoms = [ a for a in self.selected_atoms if self.target_object.atoms.get ( a.atom_id ) is a ]
        selected_ids = { a.atom_id for a in self.selected_atoms }
        try:
            self.tables = D.build_tables ( self.target_object )
        except Exception as error:
            traceback.print_exc ( )
            self.tables = None
            self.label_status.set_text ( "DYFF analysis failed: {}".format ( error ) )
            return
        for store in self.stores.values ( ): store.clear ( )
        for rid, row in enumerate ( self.tables["atoms"] ):
            values = [ ]
            for title, field, editable in _ATOM_COLUMNS:
                v = row[field]
                if field == "charge": values.append ( "" if v is None else "{:+.4f}".format ( v ) )
                elif field == "formal": values.append ( "{:+d}".format ( v ) if v else "0" )
                elif field == "atom_id": values.append ( str ( v ) )
                else: values.append ( str ( v ) )
            self.stores["atoms"].append ( [ rid ] + values + [ row["quality"], D.LEVEL_COLORS[row["level"]], row["level"],
                                                               row["atom_id"] in selected_ids ] )
        for term, title, columns in _TERM_TABS:
            for rid, row in enumerate ( self.tables[term] ):
                values = [ ]
                for name, field in columns:
                    if field == "measured":
                        m = row.get ( "measured" )
                        values.append ( "" if m is None else ( "{:.3f} A".format ( m ) if term == "bond" else "{:.1f} deg".format ( m ) ) )
                    else:
                        values.append ( str ( row.get ( field, "" ) ) )
                involved = any ( i in selected_ids for inst in row["instances"] for i in inst )
                self.stores[term].append ( [ rid ] + values + [ row["quality"], D.LEVEL_COLORS[row["level"]], row["level"], involved ] )
        self._update_tab_titles ( )
        name = self.target_object.name
        header = "{}: {} atom(s){}".format ( name, len ( self.target_object.atoms ),
                                           ", {} selected".format ( len ( selected_ids ) ) if selected_ids else "" )
        self.header_label.set_text ( header )
        partial = "{:+.4f}".format ( self.tables["total_partial"] ) if self.tables["has_charges"] else "none (DYFF: 0)"
        self.label_charge.set_markup ( "Formal charge: <b>{:+d}</b>   Partial charges: <b>{}</b>".format ( self.tables["total_formal"], partial ) )
        if self.tables["error"]:
            self.label_status.set_text ( "DYFF model not built: {} -- fix the red (untyped) atoms.".format ( self.tables["error"][:160] ) )
        self._update_labeled_atoms ( )

    def _update_labeled_atoms ( self, page = None ):
        """ Numbered labels in the 3D view for the atoms listed in the Atoms
            tab (as the Atom Types window did) while that tab is shown. """
        if not self.visible or self.tables is None: return
        page = self.notebook.get_current_page ( ) if page is None else page
        atoms = [ ]
        if page == 0:
            for row in self.filters["atoms"]:
                atom = self.target_object.atoms.get ( self.tables["atoms"][row[0]]["atom_id"] )
                if atom is not None: atoms.append ( atom )
        self.vm_session.builder_atom_types_labeled_atoms = atoms[:300]
        if getattr ( self.vm_session, "vm_glcore", None ) is not None:
            self.vm_session.vm_glcore.queue_draw ( )

    #-------------------------------------------------------------------------
    #  P R O J E C T I O N
    #-------------------------------------------------------------------------
    def on_row_selected ( self, selection, term ):
        if self.tables is None: return
        if term == "atoms":
            model, paths = selection.get_selected_rows ( )
            if not paths: return
            rid = model[paths[-1]][0]
        else:
            model, tree_iter = selection.get_selected ( )
            if tree_iter is None: return
            rid = model[tree_iter][0]
        # . only when the row really changes (a second click on the selected row
        #   starts the cell editor -- re-projecting would cancel it), and after
        #   the click is handled: centring the camera is a blocking animation
        if self.current and self.current[0] == term and self.current[1] == rid:
            return
        self.current = [ term, rid, 0 ]
        self.button_fit_torsion.set_sensitive ( term == "dihedral" )
        GLib.idle_add ( self._project_idle )

    def _project_idle ( self ):
        try:
            self._project ( )
        except Exception:
            traceback.print_exc ( )
        return False

    def on_row_activated ( self, tree, path, column ):
        """ Double-click: edit the cell under the pointer when it is editable. """
        for renderer in column.get_cells ( ):
            if isinstance ( renderer, Gtk.CellRendererText ) and renderer.get_property ( "editable" ):
                tree.set_cursor ( path, column, True )
                return

    def on_button_fit_torsion_clicked ( self, widget ):
        """ xTB scan of the dihedral shown (row + previous/next term) and fit
            of its DYFF parameter (torsion_fit_window.py). """
        if not self.current or self.current[0] != "dihedral": return
        from gui.windows.builder.torsion_fit_window import TorsionFitWindow
        term, rid, k = self.current
        row = self.tables[term][rid]
        if getattr ( self, "torsion_fit_window", None ) is None:
            self.torsion_fit_window = TorsionFitWindow ( main = self.main, on_applied = self._on_torsion_fit_applied )
        self.torsion_fit_window.open_window ( D.DYFFTorsionTarget ( self.target_object, row, row["instances"][k] ) )

    def _on_torsion_fit_applied ( self, message ):
        if self.visible:
            self.refresh ( )
            self.label_status.set_text ( message )

    def on_button_next_instance_clicked ( self, widget ): self._step ( +1 )
    def on_button_prev_instance_clicked ( self, widget ): self._step ( -1 )

    def _step ( self, delta ):
        if not self.current: return
        term, rid, k = self.current
        n = len ( self.tables[term][rid]["instances"] ) or 1
        self.current[2] = ( k + delta ) % n
        self._project ( )

    def _project ( self ):
        import numpy as np
        term, rid, k = self.current
        row = self.tables[term][rid]
        instances = row["instances"]
        if not instances: return
        atoms = [ self.target_object.atoms.get ( i ) for i in instances[k] ]
        if any ( a is None for a in atoms ): return
        text = "Term {}/{}: {}".format ( k + 1, len ( instances ), "-".join ( a.name for a in atoms ) )
        if term in ( "bond", "angle", "dihedral", "outofplane" ):
            points = [ a.coords ( ) for a in atoms ]
            measured, strain = D.term_strain ( term, row["value"], points )
            unit = "A" if term == "bond" else "deg"
            text += "   {:.3f} {}   strain {:.2f} kJ/mol".format ( measured, unit, strain ) if term == "bond" else \
                    "   {:.1f} {}   strain {:.2f} kJ/mol".format ( measured, unit, strain )
        elif term == "atoms":
            text += "   type {}   ({})".format ( row["type"], row["quality"] )
        self.label_instance.set_text ( text )
        if not self.checkbox_project.get_active ( ): return
        session = self.vm_session
        frame = getattr ( session, "selection_box_frame", None )
        if frame is not None: frame.change_toggle_button_selecting_mode_status ( True )
        else:                 session.picking_selection_mode = True
        picking = session.picking_selections
        picking.selection_function_picking ( None )
        for atom in atoms[:4]: picking.selection_function_picking ( atom )
        if self.checkbox_center.get_active ( ):
            centre = np.mean ( np.array ( [ a.coords ( ) for a in atoms ], dtype = np.float32 ), axis = 0 ).astype ( np.float32 )
            try:
                session.vm_glcore.center_on_coordinates ( self.target_object, centre )
            except Exception:
                traceback.print_exc ( )
        session.vm_glcore.queue_draw ( )

    #-------------------------------------------------------------------------
    #  A T O M   E D I T S   (applied immediately, like the Atom Types window)
    #-------------------------------------------------------------------------
    def _atom_of_path ( self, path ):
        model = self.filters["atoms"]
        child = model.convert_path_to_child_path ( Gtk.TreePath ( path ) )
        row = self.tables["atoms"][self.stores["atoms"][child][0]]
        return self.target_object.atoms.get ( row["atom_id"] ), row

    def on_atom_cell_edited ( self, renderer, path, text, field, column ):
        from gui.windows.builder.atom_types import set_manual_atom_type_override, valid_type_labels_for_symbol
        from gui.windows.builder.empty_object import sync_pdynamo_system
        atom, row = self._atom_of_path ( path )
        if atom is None: return
        text = text.strip ( )
        if field == "type":
            if not text:
                set_manual_atom_type_override ( self.target_object, atom.atom_id, None )
                message = "Atom #{}: back to the perceived type.".format ( atom.atom_id )
            else:
                valid = valid_type_labels_for_symbol ( self.vm_session, atom.symbol )
                if valid and text not in valid:
                    self.label_status.set_text ( '"{}" is not a DYFF type for {}: {}.'.format ( text, atom.symbol, ", ".join ( valid ) ) )
                    return
                set_manual_atom_type_override ( self.target_object, atom.atom_id, text )
                message = 'Atom #{}: type set to "{}".'.format ( atom.atom_id, text )
        elif field == "formal":
            try:
                atom.formal_charge = int ( text )
            except ValueError:
                self.label_status.set_text ( '"{}" is not an integer charge.'.format ( text ) ); return
            sync_pdynamo_system ( self.target_object )
            message = "Atom #{}: formal charge {:+d}.".format ( atom.atom_id, atom.formal_charge )
        else:
            try:
                value = float ( text )
            except ValueError:
                self.label_status.set_text ( '"{}" is not a number.'.format ( text ) ); return
            D.set_partial_charge ( self.target_object, atom.atom_id, value )
            message = "Atom #{}: partial charge {:+.4f}.".format ( atom.atom_id, value )
        self.refresh ( )
        self.label_status.set_text ( message )

    def on_term_cell_edited ( self, renderer, path, text, term, field ):
        """ Bond b0/k or the cosine terms "c n; c n" of a DYFF parameter (kJ/mol);
            an empty cell removes the edit (back to the DYFF rules). Applied at once
            to every term with this key (atom types + bond tags). """
        model = self.filters[term]
        child = model.convert_path_to_child_path ( Gtk.TreePath ( path ) )
        row = self.tables[term][self.stores[term][child][0]]
        text = text.strip ( )
        try:
            if not text:
                value = None
            elif term == "bond":
                b0, k = row["value"]
                if field == "v0": b0 = float ( text )
                else:             k  = float ( text )
                value = ( b0, k )
            else:
                value = D.parse_pairs ( text )
        except ValueError as error:
            self.label_status.set_text ( str ( error ) ); return
        D.set_parameter_edit ( self.target_object, term, row["key"], value )
        self.refresh ( )
        self.label_status.set_text ( "{} {}: {}.".format ( term, row["key"], "back to the DYFF rules" if value is None
                                     else "edited ({} term(s))".format ( row["count"] ) ) )

    def _selected_atom_objects ( self ):
        model, paths = self.trees["atoms"].get_selection ( ).get_selected_rows ( )
        atoms = [ ]
        for path in paths:
            atom = self.target_object.atoms.get ( self.tables["atoms"][model[path][0]]["atom_id"] )
            if atom is not None: atoms.append ( atom )
        return atoms

    def on_protonate ( self, button ):
        from gui.windows.builder.atom_ops import protonate_atom
        atoms = self._selected_atom_objects ( )
        if not atoms: self.label_status.set_text ( "Select atoms in the Atoms tab first." ); return
        messages = [ protonate_atom ( self.target_object, a.atom_id )[1] for a in atoms ]
        self.refresh ( ); self.label_status.set_text ( "  |  ".join ( messages ) )

    def on_deprotonate ( self, button ):
        from gui.windows.builder.atom_ops import deprotonate_atom
        atoms = self._selected_atom_objects ( )
        if not atoms: self.label_status.set_text ( "Select atoms in the Atoms tab first." ); return
        messages = [ deprotonate_atom ( self.target_object, a.atom_id )[1] for a in atoms ]
        self.refresh ( ); self.label_status.set_text ( "  |  ".join ( messages ) )

    def on_clear_override ( self, button ):
        from gui.windows.builder.atom_types import set_manual_atom_type_override
        atoms = self._selected_atom_objects ( )
        for a in atoms: set_manual_atom_type_override ( self.target_object, a.atom_id, None )
        self.refresh ( ); self.label_status.set_text ( "Cleared {} type override(s).".format ( len ( atoms ) ) )

    def on_rename_residue ( self, button ):
        from gui.windows.builder.atom_types_window import AtomTypesWindow
        atoms = self._selected_atom_objects ( ) or self.selected_atoms
        residues = [ ]
        for a in atoms:
            if a.residue is not None and a.residue not in residues: residues.append ( a.residue )
        if not residues: self.label_status.set_text ( "Select atoms in the Atoms tab first." ); return
        dialog = Gtk.Dialog ( title = "Rename Residue", transient_for = self.window, modal = True )
        dialog.add_buttons ( "Cancel", Gtk.ResponseType.CANCEL, "Rename", Gtk.ResponseType.OK )
        grid = Gtk.Grid ( column_spacing = 8, row_spacing = 6, margin = 10 )
        grid.attach ( Gtk.Label ( label = "Residue(s): " + ", ".join ( "{}:{} {}".format ( r.chain.name, r.name, r.index ) for r in residues[:6] ), xalign = 0 ), 0, 0, 2, 1 )
        grid.attach ( Gtk.Label ( label = "New name:", xalign = 0 ), 0, 1, 1, 1 )
        entry = Gtk.Entry ( text = residues[0].name, max_length = 4, activates_default = True )
        grid.attach ( entry, 1, 1, 1, 1 )
        grid.attach ( Gtk.Label ( label = "Number:", xalign = 0 ), 0, 2, 1, 1 )
        spin = Gtk.SpinButton.new_with_range ( -9999, 99999, 1 )
        spin.set_value ( int ( residues[0].index ) ); spin.set_sensitive ( len ( residues ) == 1 )
        grid.attach ( spin, 1, 2, 1, 1 )
        dialog.get_content_area ( ).add ( grid ); dialog.show_all ( )
        response = dialog.run ( )
        name, number = entry.get_text ( ), ( int ( spin.get_value ( ) ) if len ( residues ) == 1 else None )
        dialog.destroy ( )
        if response != Gtk.ResponseType.OK: return
        ok, message = AtomTypesWindow.apply_residue_rename ( _types.SimpleNamespace ( target_object = self.target_object ),
                                                            residues, new_name = name, new_index = number )
        self.refresh ( ); self.label_status.set_text ( message )

    def on_cm5_charges ( self, button ):
        n = len ( self.target_object.atoms )
        formal = sum ( int ( getattr ( a, "formal_charge", 0 ) or 0 ) for a in self.target_object.atoms.values ( ) )
        dialog = Gtk.Dialog ( title = "CM5 partial charges", transient_for = self.window, modal = True )
        dialog.add_buttons ( "Cancel", Gtk.ResponseType.CANCEL, "Compute", Gtk.ResponseType.OK )
        grid = Gtk.Grid ( column_spacing = 8, row_spacing = 6, margin = 10 )
        note = "GFN1-xTB CM5 charges of the whole molecule ({} atoms), averaged over equivalent atoms.".format ( n )
        if n > 300: note += "\nLarge molecule: xTB may take several minutes."
        grid.attach ( Gtk.Label ( label = note, xalign = 0 ), 0, 0, 2, 1 )
        grid.attach ( Gtk.Label ( label = "Total charge:", xalign = 0 ), 0, 1, 1, 1 )
        spin_q = Gtk.SpinButton.new_with_range ( -20, 20, 1 ); spin_q.set_value ( formal )
        grid.attach ( spin_q, 1, 1, 1, 1 )
        grid.attach ( Gtk.Label ( label = "Multiplicity:", xalign = 0 ), 0, 2, 1, 1 )
        spin_m = Gtk.SpinButton.new_with_range ( 1, 9, 1 ); spin_m.set_value ( 1 )
        grid.attach ( spin_m, 1, 2, 1, 1 )
        grid.attach ( Gtk.Label ( label = "Scale factor:", xalign = 0 ), 0, 3, 1, 1 )
        spin_s = Gtk.SpinButton.new_with_range ( 0.5, 2.0, 0.01 ); spin_s.set_digits ( 2 ); spin_s.set_value ( 1.0 )
        grid.attach ( spin_s, 1, 3, 1, 1 )
        dialog.get_content_area ( ).add ( grid ); dialog.show_all ( )
        response = dialog.run ( )
        q, m, s = int ( spin_q.get_value ( ) ), int ( spin_m.get_value ( ) ), spin_s.get_value ( )
        dialog.destroy ( )
        if response != Gtk.ResponseType.OK: return
        self.label_status.set_text ( "Running xTB..." )
        while Gtk.events_pending ( ): Gtk.main_iteration ( )
        try:
            xtb = self.vm_session.vm_config.gl_parameters.get ( 'xtb_command' ) or None
            D.compute_cm5_charges ( self.target_object, total_charge = q, multiplicity = m, scale = s, xtb = xtb )
        except Exception as error:
            traceback.print_exc ( )
            self.main.simple_dialog.error ( msg = "CM5 charges failed:\n{}".format ( error ) )
            self.label_status.set_text ( "CM5 charges failed." ); return
        self.refresh ( )
        self.label_status.set_text ( "CM5 charges set (total {:+d}, multiplicity {}); they are used by System > Add > Force Field > DYFF and Optimize.".format ( q, m ) )

    def on_clear_charges ( self, button ):
        D.clear_partial_charges ( self.target_object )
        self.refresh ( ); self.label_status.set_text ( "Partial charges removed (DYFF: all zero)." )
