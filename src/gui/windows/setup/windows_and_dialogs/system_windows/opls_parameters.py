#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: OPLS Parameters window (quality and manual editing)
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
#      Extras > OPLS > OPLS Parameters... (also opened at the end of Prepare
#      OPLS System). One tab per kind of parameter -- charges, bonds,
#      angles, dihedrals, impropers, Lennard-Jones -- each a table with the
#      values, how many terms use them, an example, the source and a
#      coloured quality column (pdynamo/opls/editor.py). Values can be
#      edited and applied to the system. Dihedrals: "Fit with xTB scan..."
#      fits the row's parameter to an xTB scan of the dihedral shown.
#

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

import os
import traceback

from gui.widgets.custom_widgets import SystemComboBox


# . tab: ( title, [ ( column title, field, editable ) ] ). Fields: "types", "v0", "v1",
#   "terms", "count", "example", "source", "residue", "name", "index", "range".
_TAB_SPECS = [
    ( "charges",    "Charges",       [ ( "Atom", "index", False ), ( "Residue", "residue", False ), ( "Name", "name", False ),
                                        ( "Type", "types", False ), ( "Charge", "v0", True ), ( "OPLS-AA class range", "range", False ),
                                        ( "Source", "source", False ) ] ),
    ( "bond",       "Bonds",         [ ( "Types", "types", False ), ( "b0 (A)", "v0", True ), ( "k (kcal/mol/A2)", "v1", True ),
                                        ( "Terms", "count", False ), ( "Example", "example", False ), ( "Source", "source", False ) ] ),
    ( "angle",      "Angles",        [ ( "Types", "types", False ), ( "theta0 (deg)", "v0", True ), ( "k (kcal/mol/rad2)", "v1", True ),
                                        ( "Terms", "count", False ), ( "Example", "example", False ), ( "Source", "source", False ) ] ),
    ( "dihedral",   "Dihedrals",     [ ( "Types", "types", False ), ( "k n phase; ...", "terms", True ),
                                        ( "Terms", "count", False ), ( "Example", "example", False ), ( "Source", "source", False ) ] ),
    ( "outofplane", "Impropers",     [ ( "Types (central first)", "types", False ), ( "K (kcal/mol)", "v0", True ),
                                        ( "Terms", "count", False ), ( "Example", "example", False ), ( "Source", "source", False ) ] ),
    ( "lj",         "Lennard-Jones", [ ( "Type", "types", False ), ( "epsilon (kcal/mol)", "v0", True ), ( "sigma (A)", "v1", True ),
                                        ( "Atoms", "count", False ), ( "Source", "source", False ) ] ),
]
_FILTERS = [ ( "ligand", "Ligand parameters" ), ( "attention", "Needs attention (yellow/red)" ), ( "all", "All parameters" ) ]


class OPLSParametersWindow:
    """ "OPLS Parameters" window. """

    def __init__ ( self, main = None ):
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        self.home       = main.home
        self.Visible    = False
        self.tables     = None
        self.pending    = { }          # term -> { key or atom index: new value }

    #-------------------------------------------------------------------------
    def open_window ( self, e_id = None ):
        if self.Visible:
            if e_id is not None: self.combobox_systems.set_active_system ( e_id = e_id )
            self.window.present ( )
            return
        self.builder = Gtk.Builder ( )
        self.builder.add_from_file ( os.path.join ( self.home, 'src/gui/windows/setup/opls_parameters_window.glade' ) )
        self.builder.connect_signals ( self )
        get = self.builder.get_object
        self.window        = get ( 'window' )
        self.notebook      = get ( 'notebook' )
        self.combo_filter  = get ( 'combo_filter' )
        self.label_charge  = get ( 'label_charge' )
        self.label_status  = get ( 'label_status' )
        self.label_legend  = get ( 'label_legend' )
        self.label_instance  = get ( 'label_instance' )
        self.checkbox_project = get ( 'checkbox_project' )
        self.checkbox_center  = get ( 'checkbox_center' )
        self.button_fit_torsion = get ( 'button_fit_torsion' )
        self.current = None                # ( term, row id, instance index )
        from pdynamo.opls.editor import LEVEL_COLORS, LEVEL_TEXT
        self.label_legend.set_markup ( "  ".join ( '<span background="{}">  </span> {}'.format ( LEVEL_COLORS[l], LEVEL_TEXT[l] )
                                                  for l in ( "green", "yellow", "red", "blue" ) ) )
        for key, text in _FILTERS:
            self.combo_filter.append ( key, text )
        self.combobox_systems = SystemComboBox ( self.main )
        get ( 'box_system' ).pack_start ( self.combobox_systems, False, False, 0 )
        self._build_tabs ( )
        self.combobox_systems.connect ( "changed", self.on_system_changed )
        self.window.connect ( 'destroy', self.close_window )
        self.window.show_all ( )
        self.Visible = True
        target = e_id if e_id is not None else self.p_session.active_id
        if self.p_session.psystem:
            self.combobox_systems.set_active_system ( e_id = target )
        self.on_system_changed ( )

    def close_window ( self, *args ):
        if not self.Visible: return
        if getattr ( self, "torsion_fit_window", None ) is not None:
            self.torsion_fit_window.close_window ( )
        self.window.destroy ( )
        self.Visible, self.tables, self.pending = False, None, { }

    def on_button_close_clicked ( self, widget ):
        self.close_window ( )

    #-------------------------------------------------------------------------
    #  T A B S
    #-------------------------------------------------------------------------
    def _build_tabs ( self ):
        """ One ListStore per tab: [ row id, <one str per column>, quality text, colour, level, ligand ]. """
        self.stores, self.filters, self.columns = { }, { }, { }
        for term, title, columns in _TAB_SPECS:
            n = len ( columns )
            store = Gtk.ListStore ( *( [ int ] + [ str ] * n + [ str, str, str, bool ] ) )
            model = store.filter_new ( )
            model.set_visible_func ( self._row_visible )
            tree = Gtk.TreeView ( model = model )
            for k, ( name, field, editable ) in enumerate ( columns ):
                renderer = Gtk.CellRendererText ( )
                if editable:
                    renderer.set_property ( "editable", True )
                    renderer.connect ( "edited", self.on_cell_edited, term, k + 1, field )
                    renderer.set_property ( "weight", 600 )
                column = Gtk.TreeViewColumn ( name, renderer, text = k + 1 )
                column.set_resizable ( True )
                column.set_sort_column_id ( k + 1 )
                tree.append_column ( column )
            renderer = Gtk.CellRendererText ( )
            renderer.set_property ( "foreground", "#1e1e1e" )        # readable on every level colour
            column = Gtk.TreeViewColumn ( "Quality", renderer, text = n + 1, background = n + 2 )
            column.set_sort_column_id ( n + 3 )
            tree.append_column ( column )
            tree.get_selection ( ).connect ( "changed", self.on_row_selected, term )
            tree.connect ( "row-activated", self.on_row_activated )
            self.trees = getattr ( self, "trees", { } )
            self.trees[term] = tree
            scrolled = Gtk.ScrolledWindow ( )
            scrolled.add ( tree )
            label = Gtk.Label ( label = title )
            self.notebook.append_page ( scrolled, label )
            self.stores[term], self.filters[term], self.columns[term] = store, model, columns
            self.tab_labels = getattr ( self, "tab_labels", { } )
            self.tab_labels[term] = ( label, title )

    def _row_visible ( self, model, tree_iter, data = None ):
        key = self.combo_filter.get_active_id ( ) or "all"
        n = model.get_n_columns ( )
        level, ligand = model[tree_iter][n - 2], model[tree_iter][n - 1]
        if key == "ligand":    return ligand
        if key == "attention": return level in ( "yellow", "red" )
        return True

    def on_combo_filter_changed ( self, widget ):
        if hasattr ( self, "filters" ):
            for model in self.filters.values ( ): model.refilter ( )
            self._update_tab_titles ( )

    def _update_tab_titles ( self ):
        """ "Bonds (shown)" followed by red/yellow counters. """
        for term, ( label, title ) in self.tab_labels.items ( ):
            store = self.stores[term]
            n = store.get_n_columns ( )
            red    = sum ( 1 for row in store if row[n - 2] == "red" )
            yellow = sum ( 1 for row in store if row[n - 2] == "yellow" )
            markup = "{} ({})".format ( title, len ( self.filters[term] ) )
            if red:    markup += ' <span foreground="#e5534b">\u25cf{}</span>'.format ( red )
            if yellow: markup += ' <span foreground="#c9a227">\u25cf{}</span>'.format ( yellow )
            label.set_markup ( markup )

    #-------------------------------------------------------------------------
    #  F I L L
    #-------------------------------------------------------------------------
    @staticmethod
    def _fmt ( value, digits = 4 ):
        return "" if value is None else ( "{:.%df}" % digits ).format ( float ( value ) )

    def _fill ( self ):
        from pdynamo.opls.editor import LEVEL_COLORS
        self.current = None
        self.button_fit_torsion.set_sensitive ( False )
        for store in self.stores.values ( ): store.clear ( )
        if self.tables is None: return
        for term, title, columns in _TAB_SPECS:
            rows = self.tables[term]
            for rid, row in enumerate ( rows ):
                values = [ ]
                for ( name, field, editable ) in columns:
                    if field == "types":
                        values.append ( row["type"] if term == "charges" else " - ".join ( row["key"] ) )
                    elif field == "v0":
                        v = row["value"]
                        if term == "charges":      values.append ( self._fmt ( v, 4 ) )
                        elif term == "outofplane": values.append ( self._fmt ( v, 3 ) )
                        else:                      values.append ( self._fmt ( v[0] if v else None, 4 ) )
                    elif field == "v1":
                        v = row["value"]
                        values.append ( self._fmt ( v[1] if v else None, 3 ) )
                    elif field == "terms":
                        v = row["value"] or [ ]
                        values.append ( "; ".join ( "{:.4f} {:d} {:.1f}".format ( k, int ( n ), p ) for ( k, n, p ) in v if float ( k ) != 0.0 ) or "none" )
                    elif field == "index":
                        values.append ( str ( row["index"] + 1 ) )
                    else:
                        values.append ( str ( row.get ( field, "" ) ) )
                self.stores[term].append ( [ rid ] + values + [ row["quality"], LEVEL_COLORS[row["level"]], row["level"], row["ligand"] ] )
        self._update_tab_titles ( )
        self._update_charge_label ( )

    def on_system_changed ( self, *args ):
        from pdynamo.opls import editor
        self.pending = { }
        e_id = self.combobox_systems.get_system_id ( )
        system = self.p_session.psystem.get ( e_id ) if e_id is not None else None
        self.system = system
        self.tables = None
        if system is None:
            self._fill ( ); return
        try:
            self.tables = editor.build_tables ( system )
            has_ligand = any ( r["ligand"] for r in self.tables["charges"] )
            if self.combo_filter.get_active_id ( ) is None:
                self.combo_filter.set_active_id ( "ligand" if has_ligand else "attention" )
            self.label_status.set_text ( "Parameter set: {}".format ( self.tables["folder"] ) )
        except Exception as error:
            self.label_status.set_text ( str ( error ) )
        self._fill ( )

    #-------------------------------------------------------------------------
    #  P R O J E C T I O N   O N   T H E   S T R U C T U R E   (picking pk1-pk4)
    #-------------------------------------------------------------------------
    def _vobject ( self ):
        """ The vismol object of the system: the visible one, else the last. """
        e_id = getattr ( self.system, "e_id", None )
        candidates = [ v for v in self.vm_session.vm_objects_dic.values ( ) if getattr ( v, "e_id", None ) == e_id ]
        if not candidates: return None
        active = [ v for v in candidates if getattr ( v, "active", False ) ]
        return ( active or candidates )[-1]

    def on_row_selected ( self, selection, term ):
        model, tree_iter = selection.get_selected ( )
        if tree_iter is None or self.tables is None: return
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
        """ xTB scan of the dihedral shown (row + previous/next term), on the
            molecule that holds it, and fit of its OPLS parameter
            (gui/windows/builder/torsion_fit_window.py, pdynamo/opls/torsion_fit.py). """
        from gui.windows.builder.torsion_fit_window import TorsionFitWindow
        from pdynamo.opls.torsion_fit import OPLSTorsionTarget
        import numpy as np
        if not self.current or self.current[0] != "dihedral" or self.system is None: return
        if self.pending:
            self.label_status.set_text ( "Apply or revert the pending changes before fitting a torsion." )
            return
        term, rid, k = self.current
        row = self.tables[term][rid]
        # . the coordinates shown (the vismol object), when it matches the system
        coordinates = None
        vobject = self._vobject ( )
        if vobject is not None and len ( vobject.atoms ) == len ( self.system.atoms ):
            frame = min ( getattr ( self.vm_session, "frame", 0 ) or 0, len ( vobject.frames ) - 1 )
            coordinates = np.array ( vobject.frames[frame], dtype = float )
        try:
            target = OPLSTorsionTarget ( self.system, row, row["instances"][k], coordinates = coordinates )
        except Exception as error:
            traceback.print_exc ( )
            self.main.simple_dialog.error ( msg = str ( error ) )
            return
        if getattr ( self, "torsion_fit_window", None ) is None:
            self.torsion_fit_window = TorsionFitWindow ( main = self.main, on_applied = self._on_torsion_fit_applied )
        self.torsion_fit_window.open_window ( target )

    def _on_torsion_fit_applied ( self, message ):
        if self.Visible:
            self.on_system_changed ( )
            self.label_status.set_text ( message )
        try:
            self.main.bottom_notebook.status_teeview_add_new_item ( message = "OPLS " + message, system = self.system )
            self.main.refresh_main_statusbar ( )
        except Exception:
            pass

    def on_button_next_instance_clicked ( self, widget ):
        self._step_instance ( +1 )

    def on_button_prev_instance_clicked ( self, widget ):
        self._step_instance ( -1 )

    def _step_instance ( self, delta ):
        if not self.current: return
        term, rid, k = self.current
        n = len ( self.tables[term][rid].get ( "instances", [ ] ) ) or 1
        self.current[2] = ( k + delta ) % n
        self._project ( )

    def _project ( self ):
        from pdynamo.opls import editor
        import numpy as np
        if not self.current or self.tables is None: return
        term, rid, k = self.current
        row = self.tables[term][rid]
        instances = row.get ( "instances", [ ] )
        if not instances:
            self.label_instance.set_text ( "" ); return
        instance = instances[k]
        vobject = self._vobject ( )
        atoms = [ vobject.atoms.get ( i ) for i in instance ] if vobject is not None else [ ]
        names = [ a.name if a is not None else "?" for a in atoms ]
        residue = ""
        if atoms and atoms[0] is not None and atoms[0].residue is not None:
            residue = "{}{} ".format ( atoms[0].residue.name, atoms[0].residue.index )
        text = "Term {}/{}: {}{}".format ( k + 1, len ( instances ), residue, "-".join ( names ) )
        if atoms and all ( a is not None for a in atoms ):
            points = [ a.coords ( ) for a in atoms ]
            value = editor.measure ( term, points )
            parameter = row["value"]
            if term == "bond" and value is not None:
                text += "   length {:.3f} A (b0 {:.3f})".format ( value, parameter[0] ) if parameter else "   length {:.3f} A".format ( value )
            elif term == "angle" and value is not None:
                text += "   angle {:.1f} deg (theta0 {:.1f})".format ( value, parameter[0] ) if parameter else "   angle {:.1f} deg".format ( value )
            elif term == "dihedral" and value is not None:
                text += "   dihedral {:.1f} deg".format ( value )
            elif term == "outofplane" and value is not None:
                text += "   central atom {:.3f} A from the plane".format ( value )
            elif term == "charges":
                text += "   charge {:+.4f}".format ( row["value"] )
            if self.checkbox_project.get_active ( ):
                self._pick ( atoms, vobject, points )
        self.label_instance.set_text ( text )

    def _pick ( self, atoms, vobject, points ):
        """ Picking mode with the term's atoms as pk1, pk2, ... (distances,
            angle and dihedral drawn by the picking machinery). """
        import numpy as np
        session = self.vm_session
        frame = getattr ( session, "selection_box_frame", None )
        if frame is not None:
            frame.change_toggle_button_selecting_mode_status ( True )
        else:
            session.picking_selection_mode = True
        picking = session.picking_selections
        picking.selection_function_picking ( None )
        for atom in atoms[:4]:
            picking.selection_function_picking ( atom )
        if self.checkbox_center.get_active ( ):
            centre = np.mean ( np.array ( points, dtype = np.float32 ), axis = 0 ).astype ( np.float32 )
            try:
                session.vm_glcore.center_on_coordinates ( vobject, centre )
            except Exception:
                traceback.print_exc ( )
        session.vm_glcore.queue_draw ( )

    #-------------------------------------------------------------------------
    #  E D I T I N G
    #-------------------------------------------------------------------------
    def on_cell_edited ( self, renderer, path, text, term, column, field ):
        from pdynamo.opls.editor import LEVEL_COLORS
        model = self.filters[term]
        child = model.convert_path_to_child_path ( Gtk.TreePath ( path ) )
        store = self.stores[term]
        row = store[child]
        source = self.tables[term][row[0]]
        try:
            if field == "terms":
                terms = [ ]
                for chunk in text.replace ( ",", " " ).split ( ";" ):
                    parts = chunk.split ( )
                    if not parts or parts[0].lower ( ) == "none": continue
                    terms.append ( ( float ( parts[0] ), int ( parts[1] ), float ( parts[2] ) ) )
                new = terms
                shown = "; ".join ( "{:.4f} {:d} {:.1f}".format ( *t ) for t in terms ) or "none"
            else:
                value = float ( text )
                shown = self._fmt ( value, 4 if field == "v0" else 3 )
        except ( ValueError, IndexError ):
            self.label_status.set_text ( 'Invalid value "{}".'.format ( text ) )
            return
        pending = self.pending.setdefault ( term, { } )
        if term == "charges":
            pending[source["index"]] = value
        elif field == "terms":
            pending[source["key"]] = new
        elif term == "outofplane":
            pending[source["key"]] = value
        else:
            current = pending.get ( source["key"], tuple ( source["value"] ) )
            current = list ( current )
            current[0 if field == "v0" else 1] = value
            pending[source["key"]] = tuple ( current )
        row[column] = shown
        n = store.get_n_columns ( )
        row[n - 4], row[n - 3] = "edited (not applied)", LEVEL_COLORS["blue"]
        self._update_charge_label ( )
        self.label_status.set_text ( "{} pending change(s): press Apply changes.".format ( sum ( len ( v ) for v in self.pending.values ( ) ) ) )

    def _pending_total_charge ( self ):
        if self.tables is None: return None
        charges = [ r["value"] for r in self.tables["charges"] ]
        for i, q in self.pending.get ( "charges", { } ).items ( ): charges[i] = q
        return sum ( charges )

    def _update_charge_label ( self ):
        total = self._pending_total_charge ( )
        if total is None:
            self.label_charge.set_text ( "" ); return
        ok = abs ( total - round ( total ) ) <= 1.0e-4
        color = "#4caf50" if ok else "#e5534b"
        self.label_charge.set_markup ( 'Total charge: <span foreground="{}"><b>{:+.4f}</b></span>'.format ( color, total ) )

    def on_button_revert_clicked ( self, widget ):
        self.on_system_changed ( )
        self.label_status.set_text ( "Changes discarded." )

    def on_button_apply_clicked ( self, widget ):
        from pdynamo.opls import editor
        if not self.pending or self.system is None:
            self.label_status.set_text ( "Nothing to apply." ); return
        try:
            before = self.system.Energy ( log = None )
        except Exception:
            before = None
        try:
            summary = editor.apply_edits ( self.system, self.pending )
        except Exception as error:
            traceback.print_exc ( )
            self.main.simple_dialog.error ( msg = str ( error ) )
            return
        try:
            after = self.system.Energy ( log = None )
        except Exception:
            after = None
        self.on_system_changed ( )
        energy = ""
        if before is not None and after is not None:
            energy = "  Energy {:.3f} -> {:.3f} kJ/mol.".format ( before, after )
        self.label_status.set_text ( "Applied: {} charge(s), {} parameter(s){}.{}".format (
                summary["charges"], summary["terms"], ", MM model rebuilt" if summary["rebuilt"] else "", energy ) )
        try:
            self.main.bottom_notebook.status_teeview_add_new_item (
                    message = "OPLS parameters edited: {} charge(s), {} parameter(s).".format ( summary["charges"], summary["terms"] ),
                    system = self.system )
            self.main.refresh_main_statusbar ( )
        except Exception:
            pass
