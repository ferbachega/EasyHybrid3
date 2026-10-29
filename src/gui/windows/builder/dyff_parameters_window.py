#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
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
import numpy as np


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
                  ( "Hybrid.", "hybridization", False ), ( "Formal q", "formal", True ),
                  ( "Partial q", "charge", True ) ]
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
            self._reorganize_layout ( )
            self.window.connect ( 'destroy', self.close_window )
            self.window.show_all ( )
            self.builder.get_object ( 'button_apply' ).hide ( )      # DYFF edits apply immediately
            self.visible = self.Visible = True
        else:
            self.window.present ( )
        self.combo_filter.set_active_id ( "selected" if self.selected_atoms else "all" )
        self.refresh ( )
        self._on_page ( self.notebook.get_current_page ( ) )

    def close_window ( self, *args ):
        if not self.visible: return
        if getattr ( self, "torsion_fit_window", None ) is not None:
            self.torsion_fit_window.close_window ( )
        self.vm_session.builder_atom_types_labeled_atoms = [ ]
        self.window.destroy ( )
        self.visible = self.Visible = False
        self.target_object, self.selected_atoms, self.tables, self.current = None, [ ], None, None
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

    #-------------------------------------------------------------------------
    #  L A Y O U T   (2026-10-24: the actions were 9 buttons in the header row,
    #  ~2100 px wide; now a sidebar of three groups, each with its current state)
    #-------------------------------------------------------------------------
    def _reorganize_layout ( self ):
        get = self.builder.get_object
        actions = get ( 'box_atom_actions' )
        actions.set_no_show_all ( True ); actions.hide ( )
        # . notebook -> paned ( sidebar | notebook )
        main_box = self.notebook.get_parent ( )
        position = main_box.child_get_property ( self.notebook, "position" )
        main_box.remove ( self.notebook )
        paned = Gtk.Paned ( orientation = Gtk.Orientation.HORIZONTAL )
        paned.pack1 ( self._build_sidebar ( ), False, False )
        paned.pack2 ( self.notebook, True, False )
        main_box.pack_start ( paned, True, True, 0 )
        main_box.reorder_child ( paned, position )
        # . the long help text goes to the bottom of the sidebar
        help_label = get ( 'label_help' )
        help_label.get_parent ( ).remove ( help_label )
        help_label.set_max_width_chars ( 34 )
        help_label.set_xalign ( 0 )
        help_label.set_markup ( "<small>" + GLib.markup_escape_text ( help_label.get_text ( ) ) + "</small>" )
        self.sidebar.pack_end ( help_label, False, False, 0 )
        # . projection bar: shorter text, Fit only on the Dihedrals tab
        self.checkbox_project.set_label ( "Show in the 3D view" )
        self.checkbox_project.set_tooltip_text ( "Show the atoms of the selected row in the structure (picking markers pk1-pk4); "
                                                 "the arrows step through the terms of the same parameter." )
        self.notebook.connect ( "switch-page", lambda nb, page, n: self._on_page ( n ) )
        # . legend and status in the bottom row
        legend = self.label_legend
        legend.get_parent ( ).remove ( legend )
        bottom = self.label_status.get_parent ( )
        bottom.pack_start ( legend, False, False, 0 )
        bottom.reorder_child ( legend, 1 )
        self.label_status.set_ellipsize ( 3 )               # . Pango.EllipsizeMode.END
        self.label_status.set_hexpand ( True )
        self.label_status.connect ( "notify::label", lambda w, *_: w.set_tooltip_text ( w.get_text ( ) ) )
        self.window.set_default_size ( 1000, 640 )
        self.window.resize ( 1000, 640 )

    def _group ( self, title ):
        frame = Gtk.Frame ( label = " {} ".format ( title ) )
        box = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 4, margin = 6 )
        frame.add ( box )
        self.sidebar.pack_start ( frame, False, False, 0 )
        return box

    def _button ( self, label, handler, tip ):
        button = Gtk.Button ( label = label )
        button.set_tooltip_text ( tip )
        button.connect ( "clicked", handler )
        return button

    def _state_label ( self, box ):
        label = Gtk.Label ( xalign = 0 )
        label.set_line_wrap ( True )
        label.set_max_width_chars ( 34 )
        box.pack_start ( label, False, False, 0 )
        return label

    def _build_sidebar ( self ):
        self.sidebar = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 8, margin_end = 6 )
        self.sidebar.set_size_request ( 250, -1 )
        # . selected atoms
        box = self._group ( "Selected atoms" )
        self.label_atoms_state = self._state_label ( box )
        grid = Gtk.Grid ( column_spacing = 4, row_spacing = 4, column_homogeneous = True )
        self.atom_buttons = [
            self._button ( "Protonate (+1)", self.on_protonate, "Add a proton to each selected atom (formal charge +1, H placed)" ),
            self._button ( "Deprotonate (-1)", self.on_deprotonate, "Remove a proton from each selected atom (formal charge -1)" ),
            self._button ( "Reset Type", self.on_clear_override, "Back to the perceived DYFF type (removes a type typed by hand)" ),
            self._button ( "Rename Residue...", self.on_rename_residue, "Rename the residue(s) of the selected atoms" ) ]
        for k, button in enumerate ( self.atom_buttons ):
            grid.attach ( button, k % 2, k // 2, 1, 1 )
        box.pack_start ( grid, False, False, 0 )
        # . partial charges
        box = self._group ( "Partial charges" )
        self.label_charges_state = self._state_label ( box )
        row = Gtk.Box ( spacing = 4, homogeneous = True )
        row.pack_start ( self._button ( "CM5 (xTB)...", self.on_cm5_charges,
                                        "Partial charges of the whole molecule: GFN1-xTB CM5, averaged over equivalent atoms" ), True, True, 0 )
        self.button_clear_charges = self._button ( "Clear", self.on_clear_charges, "Remove every partial charge (DYFF: all zero)" )
        row.pack_start ( self.button_clear_charges, True, True, 0 )
        box.pack_start ( row, False, False, 0 )
        # . xTB refinement
        box = self._group ( "xTB refinement" )
        self.label_refine_state = self._state_label ( box )
        row = Gtk.Box ( spacing = 4, homogeneous = True )
        row.pack_start ( self._button ( "Refine...", self.on_refine_xtb,
                                        "Bonds and angles from the xTB Hessian and/or torsions of the rotatable bonds (DFT "
                                        "library or relaxed xTB scans), for the molecule or the selected atoms" ), True, True, 0 )
        self.button_clear_refine = self._button ( "Clear", self.on_clear_hessian,
                                                  "Back to the DYFF rules for bonds, angles and the refined torsions" )
        row.pack_start ( self.button_clear_refine, True, True, 0 )
        box.pack_start ( row, False, False, 0 )
        return self.sidebar

    def _on_page ( self, n ):
        term = [ "atoms" ] + [ t for t, _title, _c in _TERM_TABS ]
        term = term[n] if 0 <= n < len ( term ) else None
        self.button_fit_torsion.set_visible ( term == "dihedral" )
        self._update_sidebar_state ( )

    def _update_sidebar_state ( self ):
        if not getattr ( self, "sidebar", None ) or self.target_object is None: return
        n_sel = len ( self._selected_atom_objects ( ) ) if self.tables is not None else 0
        on_atoms = self.notebook.get_current_page ( ) == 0
        if n_sel:
            self.label_atoms_state.set_markup ( "<small>{} atom(s) selected in the Atoms tab.</small>".format ( n_sel ) )
        else:
            self.label_atoms_state.set_markup ( "<small>Select rows in the Atoms tab (Ctrl or Shift + click).</small>" )
        for button in self.atom_buttons:
            button.set_sensitive ( bool ( n_sel ) and on_atoms )
        if self.tables is not None and self.tables["has_charges"]:
            self.label_charges_state.set_markup ( "<small>Set on every atom: total <b>{:+.4f}</b> e.</small>".format (
                                                  self.tables["total_partial"] ) )
            self.button_clear_charges.set_sensitive ( True )
        else:
            self.label_charges_state.set_markup ( "<small>None: DYFF uses zero charges. Compute CM5 charges or type "
                                                  "them in the Atoms tab.</small>" )
            self.button_clear_charges.set_sensitive ( False )
        edits = getattr ( self.target_object, "dyff_parameter_edits", None ) or { }
        nb, na, nt = len ( edits.get ( "bond_terms" ) or { } ), len ( edits.get ( "angle_terms" ) or { } ), \
                     len ( edits.get ( "dihedral_bonds" ) or { } )
        if nb or na or nt:
            parts = [ ]
            if nb or na: parts.append ( "Hessian: {} bond(s), {} angle(s)".format ( nb, na ) )
            if nt: parts.append ( "{} torsion(s) fitted".format ( nt ) )
            self.label_refine_state.set_markup ( "<small>{}.</small>".format ( "; ".join ( parts ) ) )
            self.button_clear_refine.set_sensitive ( True )
        else:
            self.label_refine_state.set_markup ( "<small>None: bonds, angles and torsions follow the DYFF rules.</small>" )
            self.button_clear_refine.set_sensitive ( False )

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
        self.label_charge.set_markup ( "Formal charge: <b>{:+d}</b>".format ( self.tables["total_formal"] ) )   # . partial: sidebar
        if self.tables["error"]:
            self.label_status.set_text ( "DYFF model not built: {} -- fix the red (untyped) atoms.".format ( self.tables["error"][:160] ) )
        self._update_labeled_atoms ( )
        self._update_sidebar_state ( )

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
            self._update_sidebar_state ( )
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
        if not self.current or self.tables is None: return
        term, rid, k = self.current
        if rid >= len ( self.tables.get ( term ) or [ ] ): return          # the tables changed (refresh) since the click
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
        self.label_status.set_text ( "CM5 charges set (total {:+d}, multiplicity {}); they are used by DYFF (System > Force Field) and Optimize.".format ( q, m ) )

    def _run_hessian ( self, target, q, m, solvent, adopt_geometry ):
        """ xTB optimization + Hessian -> bond/angle edits, in a worker thread. """
        import threading
        n = len ( target.atoms )
        self.label_status.set_text ( "Running xTB optimization + Hessian ({} atoms)...".format ( n ) )
        threads = max ( 1, min ( 8, ( os.cpu_count ( ) or 2 ) // 2 ) )
        def work ( ):
            try:
                result = D.compute_hessian_edits ( target, q, m, threads = threads, solvent = solvent )
                GLib.idle_add ( self._hessian_done, target, result, adopt_geometry, None )
            except Exception as error:
                traceback.print_exc ( )
                GLib.idle_add ( self._hessian_done, target, None, adopt_geometry, str ( error ) )
        threading.Thread ( target = work, daemon = True ).start ( )

    def _hessian_done ( self, target, result, adopt_geometry, error ):
        if error:
            self.main.simple_dialog.error ( msg = "xTB Hessian failed:\n{}".format ( error[-800:] ) )
            if self.visible: self.label_status.set_text ( "xTB Hessian failed." )
            return False
        try:
            summary = D.store_hessian_edits ( target, result, adopt_geometry )
        except Exception as error:
            self.main.simple_dialog.error ( msg = str ( error ) ); return False
        if getattr ( self.vm_session, "vm_glcore", None ) is not None:
            self.vm_session.vm_glcore.queue_draw ( )
        if self.visible and self.target_object is target:
            self.refresh ( )
            text = "Hessian parameters: {} bonds, {} angles (Undo reverts){}.".format ( summary["bonds"], summary["angles"],
                   "; geometry = xTB minimum" if adopt_geometry else "" )
            if summary["skipped"]:
                text += " {} term(s) kept the DYFF value (non-physical Hessian constant: {}).".format (
                        len ( summary["skipped"] ), ", ".join ( "{} {}".format ( k, "-".join ( str ( i ) for i in atoms ) )
                                                               for k, atoms, _ in summary["skipped"][:4] ) )
            if summary["n_imaginary"]:
                text += " Warning: {} imaginary frequency(ies) -- not a true minimum.".format ( summary["n_imaginary"] )
            calibrated = sum ( ( result.get ( "calibrated" ) or { } ).values ( ) )
            if calibrated:
                text += " {} term(s) calibrated to DFT by chemical environment.".format ( calibrated )
            self.label_status.set_text ( text )
        return False

    def on_refine_xtb ( self, button ):
        """ Hessian bonds/angles + per-bond torsions from relaxed xTB scans of the
            molecule (dyff_parameters.refine_torsions); worker thread. """
        import threading
        target = self.target_object
        bonds = D.rotatable_bonds ( target )
        n = len ( target.atoms )
        formal = sum ( int ( getattr ( a, "formal_charge", 0 ) or 0 ) for a in target.atoms.values ( ) )
        has_hessian = bool ( ( getattr ( target, "dyff_parameter_edits", None ) or { } ).get ( "bond_terms" ) )
        dialog = Gtk.Dialog ( title = "Refine DYFF with xTB", transient_for = self.window, modal = True )
        dialog.add_buttons ( "Cancel", Gtk.ResponseType.CANCEL, "Run", Gtk.ResponseType.OK )
        grid = Gtk.Grid ( column_spacing = 8, row_spacing = 6, margin = 10 )
        note = ( "{} atoms, {} rotatable bond(s). The Hessian gives every bond and angle its own force constant\n"
                 "(modified Seminario); each rotatable bond gets a torsion series fitted to a relaxed GFN2-xTB scan\n"
                 "(24 points) or taken from the DFT library. Selected atoms: a capped subsystem around them, whose\n"
                 "charge is the sum of its formal charges. One Undo step." ).format ( n, len ( bonds ) )
        grid.attach ( Gtk.Label ( label = note, xalign = 0 ), 0, 0, 2, 1 )
        hess = Gtk.CheckButton ( label = "Bonds and angles from the xTB Hessian (recommended)" )
        hess.set_active ( not has_hessian ); grid.attach ( hess, 0, 1, 2, 1 )
        # . 2026-10-24: the separate "Bonds/angles from xTB Hessian..." button is now these two options
        tors = Gtk.CheckButton ( label = "Torsions of the {} rotatable bond(s) (DFT library or relaxed xTB scans)".format ( len ( bonds ) ) )
        tors.set_active ( bool ( bonds ) ); grid.attach ( tors, 0, 2, 2, 1 )
        adopt = Gtk.CheckButton ( label = "Use the xTB optimized geometry (Hessian without torsions)" )
        adopt.set_active ( True ); grid.attach ( adopt, 0, 4, 2, 1 )
        def update_options ( *args ):
            if sel_radio.get_active ( ):                 # . the selection refinement always fits the torsions
                tors.set_active ( True )
            tors.set_sensitive ( not sel_radio.get_active ( ) )
            adopt.set_sensitive ( hess.get_active ( ) and not tors.get_active ( ) and not sel_radio.get_active ( ) )
            library.set_sensitive ( tors.get_active ( ) and D.torsion_library_available ( ) )
        hess.connect ( "toggled", update_options ); tors.connect ( "toggled", update_options )
        library = Gtk.CheckButton ( label = "From the DFT torsion library when the environment is known (xTB scans for the rest)" )
        library.set_tooltip_text ( "pdynamo/bond_library.py: B3LYP-D3BJ/DZVP torsion scans (QCArchive) fitted on the same DYFF + Hessian base, "
                                   "looked up by chemical environment (exact, then 2 bonds, then 1 bond around the bond). On held-out "
                                   "molecules as good as the molecule's own xTB scan (2 bonds) or better (1 bond), without scanning." )
        library.set_active ( True ); library.set_sensitive ( D.torsion_library_available ( ) )
        library.set_margin_start ( 24 ); grid.attach ( library, 0, 3, 2, 1 )
        grid.attach ( Gtk.Label ( label = "Total charge:", xalign = 0 ), 0, 8, 1, 1 )
        spin_q = Gtk.SpinButton.new_with_range ( -20, 20, 1 ); spin_q.set_value ( formal ); grid.attach ( spin_q, 1, 8, 1, 1 )
        grid.attach ( Gtk.Label ( label = "Multiplicity:", xalign = 0 ), 0, 9, 1, 1 )
        spin_m = Gtk.SpinButton.new_with_range ( 1, 9, 1 ); spin_m.set_value ( 1 ); grid.attach ( spin_m, 1, 9, 1, 1 )
        solvent = Gtk.CheckButton ( label = "Implicit water (ALPB) -- recommended for ions" ); solvent.set_active ( any ( int ( getattr ( a, "formal_charge", 0 ) or 0 ) for a in target.atoms.values ( ) ) )
        grid.attach ( solvent, 0, 10, 2, 1 )
        grid.attach ( Gtk.Label ( label = "Parallel xTB scans:", xalign = 0 ), 0, 11, 1, 1 )
        spin_w = Gtk.SpinButton.new_with_range ( 1, 32, 1 ); spin_w.set_value ( max ( 1, min ( len ( bonds ), ( os.cpu_count ( ) or 2 ) // 2 ) ) )
        grid.attach ( spin_w, 1, 11, 1, 1 )
        # . region: whole molecule or the selected atoms (capped subsystem)
        selected = self._selection_ids ( )
        whole_radio = Gtk.RadioButton.new_with_label ( None, "Whole molecule" )
        sel_radio = Gtk.RadioButton.new_with_label_from_widget ( whole_radio, "Selected atoms only ({})".format ( len ( selected ) ) )
        sel_radio.set_sensitive ( bool ( selected ) )
        grid.attach ( whole_radio, 0, 5, 1, 1 ); grid.attach ( sel_radio, 1, 5, 1, 1 )
        grid.attach ( Gtk.Label ( label = "Environment shell (bonds):", xalign = 0 ), 0, 6, 1, 1 )
        spin_shell = Gtk.SpinButton.new_with_range ( 1, 8, 1 ); spin_shell.set_value ( D.SELECTION_SHELL )
        grid.attach ( spin_shell, 1, 6, 1, 1 )
        sub_label = Gtk.Label ( label = "", xalign = 0 ); grid.attach ( sub_label, 0, 7, 2, 1 )
        def update_sub ( *args ):
            on = sel_radio.get_active ( )
            spin_shell.set_sensitive ( on )
            if not on:
                sub_label.set_text ( "" ); return
            try:
                sub = D.selection_subsystem ( target, selected, shell = int ( spin_shell.get_value ( ) ) )
                rot = sum ( 1 for b, q in bonds if b[0] in sub["selection_parent"] and b[1] in sub["selection_parent"] )
                sub_label.set_text ( "xTB subsystem: {} atoms ({} capping H, charge {:+d}); {} rotatable bond(s) inside the selection.\n"
                                     "Only bonds, angles and torsions with ALL atoms selected are refined.".format (
                                     len ( sub["atomic_numbers"] ), sub["n_caps"], sub["charge"], rot ) )
            except Exception as error:
                sub_label.set_text ( "subsystem: {}".format ( error ) )
        sel_radio.connect ( "toggled", update_sub ); spin_shell.connect ( "value-changed", update_sub )
        sel_radio.connect ( "toggled", update_options )
        if selected and len ( selected ) < n: sel_radio.set_active ( True )
        dialog.get_content_area ( ).add ( grid ); dialog.show_all ( ); update_sub ( ); update_options ( )
        response = dialog.run ( )
        options = dict ( q = int ( spin_q.get_value ( ) ), m = int ( spin_m.get_value ( ) ), hessian = hess.get_active ( ),
                         solvent = "water" if solvent.get_active ( ) else None, workers = int ( spin_w.get_value ( ) ),
                         selection = selected if sel_radio.get_active ( ) else None, shell = int ( spin_shell.get_value ( ) ),
                         library = library.get_active ( ) and library.get_sensitive ( ),
                         torsions = tors.get_active ( ), adopt = adopt.get_active ( ) )
        dialog.destroy ( )
        if response != Gtk.ResponseType.OK: return
        if not options["torsions"] and not options["selection"]:
            if not options["hessian"]:
                self.label_status.set_text ( "Nothing to refine: choose the Hessian and/or the torsions." ); return
            self._run_hessian ( target, options["q"], options["m"], options["solvent"], options["adopt"] ); return
        if options["selection"]:
            self._run_selection_refinement ( target, options ); return
        self._refine_running = True
        def say ( text ):
            GLib.idle_add ( lambda: ( self.visible and self.label_status.set_text ( "Refine with xTB: " + text ) ) and False )
        def work ( ):
            try:
                hessian_result, base = None, getattr ( target, "dyff_parameter_edits", None )
                if options["hessian"]:
                    say ( "xTB Hessian (bonds and angles)..." )
                    hessian_result = D.compute_hessian_edits ( target, options["q"], options["m"], threads = 2, solvent = options["solvent"] )
                    import copy
                    base = copy.deepcopy ( base or { } )
                    base["bond_terms"] = hessian_result["edits"]["bond_terms"]; base["angle_terms"] = hessian_result["edits"]["angle_terms"]
                    base.pop ( "term_signature", None )
                results = D.refine_torsions ( target, options["q"], options["m"], options["solvent"], workers = options["workers"],
                                              base_edits = base, progress = say, use_library = options["library"] )
                GLib.idle_add ( self._refine_done, target, results, hessian_result, None )
            except Exception as error:
                traceback.print_exc ( )
                GLib.idle_add ( self._refine_done, target, None, None, str ( error ) )
        threading.Thread ( target = work, daemon = True ).start ( )

    def _selection_ids ( self ):
        """ Atom ids selected in the Atoms tab, else the Builder selection the window was opened with. """
        ids = [ a.atom_id for a in self._selected_atom_objects ( ) ] if hasattr ( self, "trees" ) else [ ]
        if not ids: ids = [ a.atom_id for a in self.selected_atoms if self.target_object.atoms.get ( a.atom_id ) is a ]
        return sorted ( set ( ids ) )

    def _run_selection_refinement ( self, target, options ):
        import threading
        def say ( text ):
            GLib.idle_add ( lambda: ( self.visible and self.label_status.set_text ( "Refine selection with xTB: " + text ) ) and False )
        def work ( ):
            try:
                result = D.refine_selection ( target, options["selection"], multiplicity = options["m"], solvent = options["solvent"],
                                              shell = options["shell"], do_hessian = options["hessian"], workers = options["workers"],
                                              progress = say, use_library = options["library"] )
                GLib.idle_add ( self._selection_done, target, result, None )
            except Exception as error:
                traceback.print_exc ( )
                GLib.idle_add ( self._selection_done, target, None, str ( error ) )
        threading.Thread ( target = work, daemon = True ).start ( )

    def _selection_done ( self, target, result, error ):
        if error:
            self.main.simple_dialog.error ( msg = "Refine selection with xTB failed:\n{}".format ( error[-800:] ) ); return False
        self.last_refine_results = result["results"]
        counts = D.store_selection_refinement ( target, result )
        sub = result["subsystem"]
        failed = sum ( 1 for r in result["results"] if "error" in r )
        from_library = sum ( 1 for r in result["results"] if r.get ( "source" ) == "library" )
        if self.visible and self.target_object is target:
            self.refresh ( )
            self.label_status.set_text ( "Selection refined on a {}-atom subsystem ({} caps): {} bonds, {} angles, {} torsion(s){}{}{} "
                                         "(Undo reverts).".format ( sub["n_atoms"], sub["n_caps"], counts["bonds"], counts["angles"],
                                         counts["torsions"], " ({} from the DFT library)".format ( from_library ) if from_library else "",
                                         ", {} failed".format ( failed ) if failed else "",
                                         "; {} non-physical Hessian term(s) kept DYFF".format ( len ( result["skipped"] ) ) if result["skipped"] else "" ) )
        return False

    def _refine_done ( self, target, results, hessian_result, error ):
        self._refine_running = False
        if error:
            self.main.simple_dialog.error ( msg = "Refine with xTB failed:\n{}".format ( error[-800:] ) ); return False
        self.last_refine_results = results
        n = D.store_torsion_edits ( target, results, hessian_result )
        failed = [ r for r in results if "error" in r ]
        scanned = [ r for r in results if "pairs" in r and "fit_rmse" in r ]
        from_library = [ r for r in results if r.get ( "source" ) == "library" ]
        lines = [ "bond {}-{}: residual {:.2f} kcal/mol (before {:.2f}), xTB barrier {:.1f}".format (
                  r["bond"][0] + 1, r["bond"][1] + 1, r["fit_rmse"], r["before_rmse"], r["barrier_xtb"] ) for r in scanned ]
        lines += [ "bond {}-{}: DFT library ({} environment, {} scan(s))".format (
                   r["bond"][0] + 1, r["bond"][1] + 1, r["level"], r["n_scans"] ) for r in from_library ]
        self.last_refine_report = "\n".join ( lines )
        if self.visible and self.target_object is target:
            self.refresh ( )
            self.label_status.set_text ( "Refine with xTB: {}{} torsion(s){}{} (Undo reverts). {}".format (
                "bonds/angles from the Hessian + " if hessian_result else "", n,
                " ({} from the DFT library, {} from xTB scans)".format ( len ( from_library ), len ( scanned ) ) if from_library else "",
                ", {} failed".format ( len ( failed ) ) if failed else "",
                "Mean xTB-fit residual {:.2f} kcal/mol.".format ( np.mean ( [ r["fit_rmse"] for r in scanned ] ) ) if scanned else "" ) )
        return False

    def on_clear_hessian ( self, button ):
        if D.clear_hessian_edits ( self.target_object ):
            self.refresh ( ); self.label_status.set_text ( "xTB refinement removed: bonds, angles and torsions back to the DYFF rules." )
        else:
            self.label_status.set_text ( "There is no xTB refinement on this molecule." )

    def on_clear_charges ( self, button ):
        D.clear_partial_charges ( self.target_object )
        self.refresh ( ); self.label_status.set_text ( "Partial charges removed (DYFF: all zero)." )
