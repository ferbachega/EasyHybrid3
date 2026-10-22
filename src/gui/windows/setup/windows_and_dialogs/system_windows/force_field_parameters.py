#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Force Field Parameters window (quality report and manual editing)
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
#      2026-10-17, user's request: after a parametrization with xTB the
#      quality of the parameters must be shown as a report, and the
#      parameters evaluated and adjusted by hand -- as the DYFF Parameters
#      window did, but redesigned (that one had too many buttons and was
#      too wide). Replaces the OPLS Parameters window (same entry points:
#      Extras > OPLS > OPLS Parameters, end of Prepare OPLS System, and the
#      OPLS parametrization of opls_assign.py, which opens it at the end).
#
#      Layout:
#        header   system, force field / atoms / total charge, a menu (report
#                 file, parameter folder);
#        sidebar  "Report" + one entry per kind of parameter, with its
#                 red/yellow counters;
#        Report   a verdict, one card per parametrized molecule and the list
#                 of the problems -- each one opens its parameter;
#        tables   a filter (Problems / Molecule / All), the parameters with a
#                 coloured quality dot, and a detail pane for the selected
#                 one: the value measured in the structure next to the
#                 equilibrium value, editable fields, the term shown in the
#                 3D view (picking pk1-pk4) with previous/next, and, for
#                 dihedrals, the xTB scan fit;
#        footer   status, Revert / Apply / Close.
#
#      The window talks to a backend (OPLSBackend below, on
#      pdynamo/opls/editor.py); another force field (DYFF) can be added as
#      a second backend with the same methods.
#
#      2026-10-21: MMInspectBackend (pdynamo/mm_inspect.py) shows, read-only,
#      the parameters of any other MM model (AMBER, CHARMM, DYFF, OPLS not
#      prepared by EasyHybrid) straight from the pDynamo term containers.
#      Opened by System > Force Field > Force Field Parameters...
#

import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk, GLib, Pango

import os
import subprocess
import traceback

LEVEL_COLORS = { "green": "#4caf50", "yellow": "#e0b030", "red": "#e5534b", "blue": "#4a90d9", "gray": "#9e9e9e" }
FILTERS = ( ( "attention", "Problems" ), ( "ligand", "Parametrized" ), ( "all", "All" ) )


def _fmt ( value, digits = 4 ):
    return "" if value is None else ( "{:.%df}" % digits ).format ( float ( value ) )


def _terms_text ( terms ):
    return "; ".join ( "{:.4f} {:d} {:.1f}".format ( float ( k ), int ( n ), float ( p ) ) for ( k, n, p ) in ( terms or [ ] ) if float ( k ) != 0.0 ) or "none"


def _parse_terms ( text ):
    terms = [ ]
    for chunk in text.replace ( ",", " " ).split ( ";" ):
        parts = chunk.split ( )
        if not parts or parts[0].lower ( ) == "none": continue
        terms.append ( ( float ( parts[0] ), int ( parts[1] ), float ( parts[2] ) ) )
    return terms


#=============================================================================
#  O P L S   B A C K E N D
#=============================================================================
class OPLSBackend:
    """ Tables, report, measurements and editing of an OPLS system prepared by
        EasyHybrid (pdynamo/opls/editor.py). Pending edits are kept here in the
        editor.apply_edits() format until apply(). """

    force_field = "OPLS-AA"
    kinds = ( ( "charges", "Charges" ), ( "bond", "Bonds" ), ( "angle", "Angles" ), ( "dihedral", "Dihedrals" ),
              ( "outofplane", "Impropers" ), ( "lj", "Lennard-Jones" ) )

    def __init__ ( self, system ):
        from pdynamo.opls import editor
        self.system  = system
        self.tables  = editor.build_tables ( system )
        self.pending = { }
        self.quality = self._read_quality ( )

    @staticmethod
    def available ( system ):
        from pdynamo.opls import editor
        return system is not None and getattr ( system, "mmModel", None ) is not None and editor.parameter_folder ( system ) is not None

    def _read_quality ( self ):
        import json
        path = os.path.join ( os.path.dirname ( self.tables["folder"] ), "quality.json" )
        try:
            with open ( path ) as handle: return json.load ( handle )
        except ( OSError, ValueError ):
            return { }

    # . files
    def report_file ( self ):
        path = ( getattr ( self.system, "e_input_files", None ) or { } ).get ( "report" )
        if path and os.path.exists ( path ): return path
        path = os.path.join ( os.path.dirname ( self.tables["folder"] ), "prepare_opls_report.md" )
        return path if os.path.exists ( path ) else None

    def folder ( self ):
        return self.tables["folder"]

    # . rows
    def rows ( self, term ):
        return self.tables[term]

    def has_molecules ( self ):
        return bool ( self.tables["ligands"] )

    def columns ( self, term ):
        """ [ ( title, function ( row ) -> text ) ] after the quality dot. """
        if term == "charges":
            return [ ( "Atom", lambda r: "{} {}".format ( r["residue"], r["name"] ) ), ( "Type", lambda r: r["type"] ),
                     ( "Charge", lambda r: self._shown ( term, r, 0 ) ), ( "Class range", lambda r: r["range"] ) ]
        types = lambda r: " - ".join ( r["key"] )
        if term == "bond":
            return [ ( "Types", types ), ( "b0 (Å)", lambda r: self._shown ( term, r, 0 ) ),
                     ( "k (kcal/mol/Å²)", lambda r: self._shown ( term, r, 1 ) ), ( "Terms", lambda r: str ( r["count"] ) ) ]
        if term == "angle":
            return [ ( "Types", types ), ( "θ0 (°)", lambda r: self._shown ( term, r, 0 ) ),
                     ( "k (kcal/mol/rad²)", lambda r: self._shown ( term, r, 1 ) ), ( "Terms", lambda r: str ( r["count"] ) ) ]
        if term == "dihedral":
            return [ ( "Types", types ), ( "k n phase", lambda r: self._shown ( term, r, 0 ) ), ( "Terms", lambda r: str ( r["count"] ) ) ]
        if term == "outofplane":
            return [ ( "Types (central first)", types ), ( "K (kcal/mol)", lambda r: self._shown ( term, r, 0 ) ),
                     ( "Terms", lambda r: str ( r["count"] ) ) ]
        return [ ( "Type", lambda r: r["key"][0] ), ( "ε (kcal/mol)", lambda r: self._shown ( term, r, 0 ) ),
                 ( "σ (Å)", lambda r: self._shown ( term, r, 1 ) ), ( "Atoms", lambda r: str ( r["count"] ) ) ]

    def _value ( self, term, row ):
        """ Current value including a pending edit. """
        pending = self.pending.get ( term, { } )
        key = row["index"] if term == "charges" else row["key"]
        return pending.get ( key, row["value"] )

    def is_pending ( self, term, row ):
        return ( row["index"] if term == "charges" else row["key"] ) in self.pending.get ( term, { } )

    def _shown ( self, term, row, slot ):
        v = self._value ( term, row )
        if term == "charges":    return _fmt ( v, 4 )
        if term == "dihedral":   return _terms_text ( v )
        if term == "outofplane": return _fmt ( v, 3 )
        return _fmt ( v[slot] if v else None, 4 if slot == 0 else 3 )

    def level ( self, term, row ):
        return "blue" if self.is_pending ( term, row ) else row["level"]

    def quality_text ( self, term, row ):
        return "edited (not applied)" if self.is_pending ( term, row ) else row["quality"]

    def source_text ( self, term, row ):
        return row.get ( "source", "" )

    def in_molecule ( self, term, row ):
        return row["ligand"]

    # . detail pane
    def title ( self, term, row ):
        if term == "charges": return "Charge of {} {} ({})".format ( row["residue"], row["name"], row["type"] )
        names = { "bond": "Bond", "angle": "Angle", "dihedral": "Dihedral", "outofplane": "Improper", "lj": "Lennard-Jones" }
        return "{} {}".format ( names[term], " - ".join ( row["key"] ) )

    def fields ( self, term, row ):
        """ Editable fields: [ ( slot, label, text, unit ) ]. """
        v = self._value ( term, row )
        if term == "charges":    return [ ( 0, "Charge", _fmt ( v, 4 ), "e" ) ]
        if term == "bond":       return [ ( 0, "b0", _fmt ( v[0], 4 ), "Å" ), ( 1, "k", _fmt ( v[1], 2 ), "kcal/mol/Å²" ) ]
        if term == "angle":      return [ ( 0, "θ0", _fmt ( v[0], 2 ), "°" ), ( 1, "k", _fmt ( v[1], 2 ), "kcal/mol/rad²" ) ]
        if term == "dihedral":   return [ ( 0, "Terms", _terms_text ( v ), "k (kcal/mol) n phase; ..." ) ]
        if term == "outofplane": return [ ( 0, "K", _fmt ( v, 3 ), "kcal/mol" ) ]
        return [ ( 0, "ε", _fmt ( v[0], 4 ), "kcal/mol" ), ( 1, "σ", _fmt ( v[1], 4 ), "Å" ) ]

    def set_field ( self, term, row, slot, text ):
        """ Stores a pending edit; raises ValueError for a bad value. Returns True if it changed. """
        pending = self.pending.setdefault ( term, { } )
        if term == "dihedral":
            new = _parse_terms ( text )
            if _terms_text ( new ) == _terms_text ( self._value ( term, row ) ): return False
            pending[row["key"]] = new
            return True
        value = float ( text )
        if term == "charges":
            if abs ( value - self._value ( term, row ) ) < 1e-6: return False
            pending[row["index"]] = value
        elif term == "outofplane":
            if abs ( value - float ( self._value ( term, row ) ) ) < 1e-6: return False
            pending[row["key"]] = value
        else:
            current = list ( self._value ( term, row ) )
            if abs ( value - float ( current[slot] ) ) < 1e-6: return False
            current[slot] = value
            pending[row["key"]] = tuple ( current )
        return True

    def n_pending ( self ):
        return sum ( len ( v ) for v in self.pending.values ( ) )

    def total_charge ( self ):
        charges = [ r["value"] for r in self.tables["charges"] ]
        for i, q in self.pending.get ( "charges", { } ).items ( ): charges[i] = q
        return sum ( charges )

    def apply ( self ):
        from pdynamo.opls import editor
        summary = editor.apply_edits ( self.system, self.pending )
        return "{} charge(s), {} parameter(s){}".format ( summary["charges"], summary["terms"],
                                                         ", MM model rebuilt" if summary["rebuilt"] else "" )

    def measure_text ( self, term, row, points ):
        from pdynamo.opls import editor
        value = editor.measure ( term, points )
        v = self._value ( term, row )
        if term == "bond" and value is not None:
            return "measured {:.3f} Å   ·   b0 {:.3f} Å".format ( value, v[0] ) if v else "measured {:.3f} Å".format ( value )
        if term == "angle" and value is not None:
            return "measured {:.1f}°   ·   θ0 {:.1f}°".format ( value, v[0] ) if v else "measured {:.1f}°".format ( value )
        if term == "dihedral" and value is not None:
            return "measured {:.1f}°".format ( value )
        if term == "outofplane" and value is not None:
            return "central atom {:.3f} Å from the plane".format ( value )
        if term == "charges":
            return "charge {:+.4f} e   ·   OPLS-AA class range {}".format ( v, row["range"] or "unknown" )
        return ""

    def can_fit_torsion ( self, term ):
        return term == "dihedral"

    def torsion_target ( self, row, instance, coordinates ):
        from pdynamo.opls.torsion_fit import OPLSTorsionTarget
        return OPLSTorsionTarget ( self.system, row, instance, coordinates = coordinates )

    # . report
    def report ( self ):
        """ { "verdict": ( level, title, text ), "molecules": [ ( name, levels, recommendations ) ],
              "issues": [ ( level, text, term, row index or None ) ] }. """
        issues = [ ]
        for term, title in self.kinds:
            for rid, row in enumerate ( self.tables[term] ):
                if row["level"] == "red":
                    issues.append ( ( "red", "{}: {} -- {}".format ( title.rstrip ( "s" ) if term != "lj" else title,
                                     self._row_name ( term, row ), row["quality"] ), term, rid ) )
        yellow = { }
        for term, title in self.kinds:
            for rid, row in enumerate ( self.tables[term] ):
                if row["level"] == "yellow":
                    yellow.setdefault ( ( term, row["quality"] ), [ ] ).append ( rid )
        for ( term, quality ), rids in yellow.items ( ):
            title = dict ( self.kinds )[term]
            if len ( rids ) == 1:
                row = self.tables[term][rids[0]]
                issues.append ( ( "yellow", "{}: {} -- {}".format ( title, self._row_name ( term, row ), quality ), term, rids[0] ) )
            else:
                issues.append ( ( "yellow", "{}: {} parameters -- {}".format ( title, len ( rids ), quality ), term, None ) )
        n_red = sum ( 1 for i in issues if i[0] == "red" )
        n_yellow = sum ( 1 if i[3] is not None else 0 for i in issues if i[0] == "yellow" ) + \
                   sum ( len ( r ) for ( t, q ), r in yellow.items ( ) if len ( r ) > 1 )
        if n_red:
            verdict = ( "red", "Needs attention", "{} parameter(s) were estimated or are missing (red). Check them before "
                        "production simulations; the recommendations below say what to do.".format ( n_red ) )
        elif n_yellow:
            verdict = ( "yellow", "Usable -- review the generic terms", "Every bond and angle comes from the OPLS database; "
                        "{} parameter(s) are generic or charges outside their class range (yellow).".format ( n_yellow ) )
        else:
            verdict = ( "green", "Ready", "Every parameter comes from the OPLS-AA database and every charge is within the "
                        "range of its class." )
        molecules = [ ( path, data.get ( "levels", { } ), data.get ( "recommendations", [ ] ) ) for path, data in self.quality.items ( ) ]
        return { "verdict": verdict, "molecules": molecules, "issues": issues }

    def _row_name ( self, term, row ):
        if term == "charges": return "{} {}".format ( row["residue"], row["name"] )
        return "{} (e.g. {})".format ( "-".join ( row["key"] ), row["example"] ) if row.get ( "example" ) else "-".join ( row["key"] )

    def summary ( self ):
        return "{} · {} atoms · {} parameter set".format ( self.force_field, len ( self.system.atoms ),
                                                        "custom" if self.has_molecules ( ) else "standard" )


#=============================================================================
#  A N Y   M M   M O D E L   ( R E A D - O N L Y )
#=============================================================================
class MMInspectBackend:
    """ 2026-10-21: inspection of the MM parameters of any system (OPLS, AMBER,
        CHARMM, DYFF...) read from its pDynamo term containers
        (pdynamo/mm_inspect.py). Nothing is editable here. """

    editable = False

    def __init__ ( self, system ):
        from pdynamo import mm_inspect
        self.system      = system
        self.tables      = mm_inspect.build_tables ( system )
        self.force_field = self.tables["force_field"]
        self.kinds       = tuple ( self.tables["kinds"] )
        self.pending     = { }

    @staticmethod
    def available ( system ):
        return system is not None and getattr ( system, "mmModel", None ) is not None and getattr ( system, "mmState", None ) is not None

    def report_file ( self ): return None
    def folder ( self ):
        files = getattr ( self.system, "e_input_files", None ) or { }
        for key in ( "opls_parameters", "amber_prmtop", "charmm_psf", "charmm_par" ):
            path = files.get ( key )
            if isinstance ( path, ( list, tuple ) ): path = path[0] if path else None
            if path and os.path.exists ( path ): return path if os.path.isdir ( path ) else os.path.dirname ( path )
        return None

    def rows ( self, term ): return self.tables[term]
    def has_molecules ( self ): return False
    def columns ( self, term ): return self.tables["columns"][term]
    def is_pending ( self, term, row ): return False
    def level ( self, term, row ): return "gray"
    def quality_text ( self, term, row ): return ""
    def source_text ( self, term, row ): return row.get ( "source", "" )
    def in_molecule ( self, term, row ): return False

    def title ( self, term, row ):
        if term == "charges": return "Charge of {} {} ({})".format ( row["residue"], row["name"], row["type"] )
        return "{} {}".format ( dict ( self.kinds )[term], " - ".join ( row["key"] ) )

    def fields ( self, term, row ):
        """ The row's values, shown read-only. """
        return [ ( k, title, func ( row ), "" ) for k, ( title, func ) in enumerate ( self.columns ( term ) )
                 if title not in ( "Types", "Atom", "Terms", "Atoms", "LJ type" ) ]

    def set_field ( self, term, row, slot, text ): return False
    def n_pending ( self ): return 0
    def total_charge ( self ): return self.tables["total_charge"]
    def apply ( self ): return ""

    def measure_text ( self, term, row, points ):
        from pdynamo import mm_inspect
        kind = self.tables["measure"].get ( term )
        value = mm_inspect.measure ( kind, points ) if kind else None
        if value is None:
            return "charge {:+.4f} e".format ( row["value"] ) if term == "charges" else ""
        if kind == "bond":       return "measured {:.3f} Å".format ( value )
        if kind == "outofplane": return "central atom {:.3f} Å from the plane".format ( value )
        return "measured {:.1f}°".format ( value )

    def can_fit_torsion ( self, term ): return False

    def report ( self ):
        total = self.tables["total_charge"]
        if abs ( total ) < 5.0e-5: total = 0.0
        lines = [ "{} atoms, total charge {:+.4f} e.".format ( len ( self.system.atoms ), total ) ]
        for term, title in self.kinds[1:]:
            rows = self.tables[term]
            unit = "types" if term == "lj" else "parameters"
            lines.append ( "{}: {} {}, {} terms.".format ( title, len ( rows ), unit, sum ( r["count"] for r in rows ) ) if term != "lj"
                           else "{}: {} types.".format ( title, len ( rows ) ) )
        lines.append ( "Units: kJ/mol and Å (pDynamo); harmonic terms are E = k (x - x0)². Read-only: select a parameter "
                       "to see where it is used in the structure." )
        issues = [ ]
        if abs ( total - round ( total ) ) > 1.0e-3:
            issues.append ( ( "yellow", "Total charge {:+.4f} e is not an integer".format ( total ), "charges", None ) )
        return { "verdict": ( "blue", "{} parameters".format ( self.force_field ), "\n".join ( lines ) ),
                 "molecules": [ ], "issues": issues }

    def summary ( self ):
        return "{} · {} atoms · read-only".format ( self.force_field, len ( self.system.atoms ) )


#=============================================================================
#  W I N D O W
#=============================================================================
class ForceFieldParametersWindow:
    """ "Force Field Parameters" window. open_window ( e_id ) as the old OPLS
        Parameters window. """

    def __init__ ( self, main = None ):
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        self.home       = main.home
        self.Visible    = False
        self.backend    = None
        self.current    = None             # [ term, row index, instance index ]

    #-------------------------------------------------------------------------
    def open_window ( self, e_id = None ):
        if self.Visible:
            if e_id is not None: self.combobox_systems.set_active_system ( e_id = e_id )
            self.window.present ( )
            return
        from gui.widgets.custom_widgets import SystemComboBox
        self.window = Gtk.Window ( title = "Force Field Parameters" )
        self.window.set_default_size ( 900, 600 )
        self.window.set_transient_for ( self.main.window )
        self.window.connect ( 'destroy', self.close_window )
        outer = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 6, margin = 8 )
        self.window.add ( outer )

        # . header
        header = Gtk.Box ( spacing = 8 )
        self.combobox_systems = SystemComboBox ( self.main )
        header.pack_start ( self.combobox_systems, False, False, 0 )
        self.label_summary = Gtk.Label ( xalign = 0 )
        self.label_summary.set_ellipsize ( Pango.EllipsizeMode.END )
        header.pack_start ( self.label_summary, True, True, 0 )
        menu_button = Gtk.MenuButton ( )
        menu_button.set_image ( Gtk.Image.new_from_icon_name ( "open-menu-symbolic", Gtk.IconSize.BUTTON ) )
        menu_button.set_tooltip_text ( "Report file and parameter folder" )
        menu = Gtk.Menu ( )
        for text, handler in ( ( "Open the report file", self.on_open_report ), ( "Open the parameter folder", self.on_open_folder ) ):
            item = Gtk.MenuItem ( label = text ); item.connect ( "activate", handler ); menu.append ( item )
        menu.show_all ( )
        menu_button.set_popup ( menu )
        header.pack_end ( menu_button, False, False, 0 )
        outer.pack_start ( header, False, False, 0 )

        # . sidebar + pages
        body = Gtk.Box ( spacing = 0 )
        self.sidebar = Gtk.ListBox ( )
        self.sidebar.set_selection_mode ( Gtk.SelectionMode.SINGLE )
        self.sidebar.connect ( "row-selected", self.on_sidebar_selected )
        side_frame = Gtk.Frame ( )
        side_frame.add ( self.sidebar )
        side_frame.set_size_request ( 175, -1 )
        body.pack_start ( side_frame, False, False, 0 )
        right = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 6, margin_start = 8 )
        body.pack_start ( right, True, True, 0 )
        outer.pack_start ( body, True, True, 0 )

        # . filter (table pages only)
        self.filter_box = Gtk.Box ( spacing = 0 )
        self.filter_box.get_style_context ( ).add_class ( "linked" )
        self.filter_buttons = { }
        group = None
        for key, text in FILTERS:
            button = Gtk.RadioButton.new_with_label_from_widget ( group, text )
            button.set_mode ( False )
            button.connect ( "toggled", self.on_filter_toggled, key )
            self.filter_box.pack_start ( button, False, False, 0 )
            self.filter_buttons[key] = button
            group = group or button
        self.label_filter_info = Gtk.Label ( xalign = 1 )
        self.label_filter_info.get_style_context ( ).add_class ( "dim-label" )
        filter_row = Gtk.Box ( spacing = 8 )
        filter_row.pack_start ( self.filter_box, False, False, 0 )
        filter_row.pack_end ( self.label_filter_info, False, False, 0 )
        self.filter_row = filter_row
        right.pack_start ( filter_row, False, False, 0 )

        self.stack = Gtk.Stack ( )
        right.pack_start ( self.stack, True, True, 0 )
        self.report_box = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 10, margin = 4 )
        scrolled = Gtk.ScrolledWindow ( )
        scrolled.add ( self.report_box )
        self.stack.add_named ( scrolled, "report" )
        self.trees, self.stores, self.filtered = { }, { }, { }

        # . detail pane
        self.detail = self._build_detail ( )
        right.pack_start ( self.detail, False, False, 0 )

        # . footer
        footer = Gtk.Box ( spacing = 6 )
        self.label_status = Gtk.Label ( xalign = 0 )
        self.label_status.set_ellipsize ( Pango.EllipsizeMode.END )
        footer.pack_start ( self.label_status, True, True, 0 )
        self.button_revert = Gtk.Button ( label = "Revert" )
        self.button_revert.set_tooltip_text ( "Discard the changes not applied yet" )
        self.button_revert.connect ( "clicked", self.on_revert )
        self.button_apply = Gtk.Button ( label = "Apply" )
        self.button_apply.get_style_context ( ).add_class ( "suggested-action" )
        self.button_apply.connect ( "clicked", self.on_apply )
        close = Gtk.Button ( label = "Close" )
        close.connect ( "clicked", lambda b: self.close_window ( ) )
        for b in ( self.button_revert, self.button_apply, close ): footer.pack_start ( b, False, False, 0 )
        outer.pack_start ( footer, False, False, 0 )

        self.combobox_systems.connect ( "changed", self.on_system_changed )
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
        self.Visible = False
        self.window.destroy ( )
        self.backend, self.current = None, None

    #-------------------------------------------------------------------------
    def _build_detail ( self ):
        frame = Gtk.Frame ( )
        box = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 6, margin = 8 )
        frame.add ( box )
        top = Gtk.Box ( spacing = 8 )
        self.detail_dot = Gtk.Label ( )
        self.detail_title = Gtk.Label ( xalign = 0 )
        self.detail_title.set_ellipsize ( Pango.EllipsizeMode.END )
        self.detail_quality = Gtk.Label ( xalign = 1 )
        self.detail_quality.get_style_context ( ).add_class ( "dim-label" )
        top.pack_start ( self.detail_dot, False, False, 0 )
        top.pack_start ( self.detail_title, True, True, 0 )
        top.pack_end ( self.detail_quality, False, False, 0 )
        box.pack_start ( top, False, False, 0 )
        self.detail_fields = Gtk.Box ( spacing = 8 )
        box.pack_start ( self.detail_fields, False, False, 0 )
        bottom = Gtk.Box ( spacing = 6 )
        self.button_prev = Gtk.Button.new_from_icon_name ( "go-previous-symbolic", Gtk.IconSize.BUTTON )
        self.button_next = Gtk.Button.new_from_icon_name ( "go-next-symbolic", Gtk.IconSize.BUTTON )
        self.button_prev.set_tooltip_text ( "Previous term with this parameter" )
        self.button_next.set_tooltip_text ( "Next term with this parameter" )
        self.button_prev.connect ( "clicked", lambda b: self._step ( -1 ) )
        self.button_next.connect ( "clicked", lambda b: self._step ( +1 ) )
        self.label_instance = Gtk.Label ( xalign = 0 )
        self.label_instance.set_ellipsize ( Pango.EllipsizeMode.END )
        self.check_show = Gtk.CheckButton ( label = "Show in 3D" )
        self.check_show.set_active ( True )
        self.check_show.set_tooltip_text ( "Mark the term's atoms in the structure (picking pk1-pk4) and centre the view on them" )
        self.check_show.connect ( "toggled", lambda b: self._project ( ) )
        self.button_fit = Gtk.Button ( label = "Fit with xTB scan..." )
        self.button_fit.set_tooltip_text ( "Scan this dihedral with xTB and fit the parameter to the scan" )
        self.button_fit.connect ( "clicked", self.on_fit_torsion )
        bottom.pack_start ( self.button_prev, False, False, 0 )
        bottom.pack_start ( self.button_next, False, False, 0 )
        bottom.pack_start ( self.label_instance, True, True, 0 )
        bottom.pack_end ( self.button_fit, False, False, 0 )
        bottom.pack_end ( self.check_show, False, False, 0 )
        box.pack_start ( bottom, False, False, 0 )
        return frame

    def _build_table ( self, term ):
        """ ListStore: [ row index, dot markup, <column texts...>, quality, level, in molecule ]. """
        columns = self.backend.columns ( term )
        n = len ( columns )
        store = Gtk.ListStore ( *( [ int, str ] + [ str ] * n + [ str, str, bool ] ) )
        model = store.filter_new ( )
        model.set_visible_func ( self._row_visible )
        tree = Gtk.TreeView ( model = model )
        tree.set_enable_search ( False )
        renderer = Gtk.CellRendererText ( )
        dots = Gtk.TreeViewColumn ( "", renderer, markup = 1 )
        dots.set_visible ( getattr ( self.backend, "editable", True ) )
        tree.append_column ( dots )
        for k, ( title, func ) in enumerate ( columns ):
            renderer = Gtk.CellRendererText ( )
            column = Gtk.TreeViewColumn ( title, renderer, text = k + 2 )
            column.set_resizable ( True )
            column.set_sort_column_id ( k + 2 )
            tree.append_column ( column )
        renderer = Gtk.CellRendererText ( )
        renderer.set_property ( "foreground", "#888888" )
        column = Gtk.TreeViewColumn ( "Quality", renderer, text = n + 2 )
        column.set_sort_column_id ( n + 3 )
        column.set_visible ( getattr ( self.backend, "editable", True ) )
        tree.append_column ( column )
        tree.get_selection ( ).connect ( "changed", self.on_row_selected, term )
        scrolled = Gtk.ScrolledWindow ( )
        scrolled.add ( tree )
        scrolled.show_all ( )
        self.stack.add_named ( scrolled, term )
        self.trees[term], self.stores[term], self.filtered[term] = tree, store, model

    #-------------------------------------------------------------------------
    #  S Y S T E M   A N D   F I L L
    #-------------------------------------------------------------------------
    def on_system_changed ( self, *args ):
        e_id = self.combobox_systems.get_system_id ( )
        system = self.p_session.psystem.get ( e_id ) if e_id is not None else None
        self.backend, self.current = None, None
        for term in list ( self.trees ):
            self.stack.remove ( self.stack.get_child_by_name ( term ) )
        self.trees, self.stores, self.filtered = { }, { }, { }
        for row in self.sidebar.get_children ( ): self.sidebar.remove ( row )
        message = None
        try:
            if system is None:
                message = "No system."
            elif OPLSBackend.available ( system ):
                self.backend = OPLSBackend ( system )
            elif MMInspectBackend.available ( system ):
                self.backend = MMInspectBackend ( system )
            else:
                message = ( "'{}' has no force field. Define one with System > Force Field > Define Force Field "
                            "(Prepare OPLS System or Prepare DYFF System).".format ( system.label ) )
        except Exception as error:
            traceback.print_exc ( )
            message = str ( error )
        if self.backend is None:
            self.label_summary.set_text ( "" )
            self._show_message ( message )
            self.stack.set_visible_child_name ( "report" )
            self.filter_row.hide ( ); self.detail.hide ( )
            self._update_footer ( )
            return
        for term, title in self.backend.kinds:
            self._build_table ( term )
        self._add_sidebar_row ( "report", "Report" )
        for term, title in self.backend.kinds:
            self._add_sidebar_row ( term, title )
        self.sidebar.show_all ( )
        editable = getattr ( self.backend, "editable", True )
        self.filter_buttons["ligand"].set_sensitive ( self.backend.has_molecules ( ) )
        self.filter_buttons["attention"].set_sensitive ( editable )
        if not editable: self.filter_buttons["all"].set_active ( True )
        self.filter_box.set_visible ( editable )
        self.button_apply.set_visible ( editable ); self.button_revert.set_visible ( editable )
        self._fill ( )
        self.sidebar.select_row ( self.sidebar.get_row_at_index ( 0 ) )

    def _add_sidebar_row ( self, name, title ):
        row = Gtk.ListBoxRow ( )
        row.name = name
        box = Gtk.Box ( spacing = 6, margin = 6 )
        label = Gtk.Label ( label = title, xalign = 0 )
        badge = Gtk.Label ( xalign = 1 )
        box.pack_start ( label, True, True, 0 )
        box.pack_end ( badge, False, False, 0 )
        row.add ( box )
        row.badge = badge
        self.sidebar.add ( row )

    def _fill ( self ):
        b = self.backend
        for term, title in b.kinds:
            store = self.stores[term]
            store.clear ( )
            columns = b.columns ( term )
            for rid, row in enumerate ( b.rows ( term ) ):
                level = b.level ( term, row )
                store.append ( [ rid, '<span foreground="{}">●</span>'.format ( LEVEL_COLORS[level] ) ]
                               + [ func ( row ) for title_, func in columns ] + [ b.quality_text ( term, row ), level, b.in_molecule ( term, row ) ] )
        self._update_badges ( )
        self._fill_report ( )
        self.label_summary.set_markup ( "<b>{}</b>".format ( GLib.markup_escape_text ( b.summary ( ) ) ) )
        if not getattr ( b, "editable", True ):
            self.filter_buttons["all"].set_active ( True )
        elif not any ( button.get_active ( ) and button.get_sensitive ( ) for button in self.filter_buttons.values ( ) ) or \
           not self.filter_buttons["ligand"].get_sensitive ( ) and self.filter_buttons["ligand"].get_active ( ):
            self.filter_buttons["attention"].set_active ( True )
        self._update_footer ( )

    def _update_badges ( self ):
        b = self.backend
        for row in self.sidebar.get_children ( ):
            if row.name == "report":
                level = b.report ( )["verdict"][0]
                row.badge.set_markup ( '<span foreground="{}">●</span>'.format ( LEVEL_COLORS[level] ) )
                continue
            levels = [ b.level ( row.name, r ) for r in b.rows ( row.name ) ]
            parts = [ ]
            for level in ( "red", "yellow", "blue" ):
                n = levels.count ( level )
                if n: parts.append ( '<span foreground="{}">●{}</span>'.format ( LEVEL_COLORS[level], n ) )
            parts.append ( '<span foreground="#888888">{}</span>'.format ( len ( levels ) ) )
            row.badge.set_markup ( " ".join ( parts ) )

    def _show_message ( self, text ):
        for child in self.report_box.get_children ( ): self.report_box.remove ( child )
        label = Gtk.Label ( label = text or "", xalign = 0, wrap = True )
        self.report_box.pack_start ( label, False, False, 0 )
        self.report_box.show_all ( )

    def _fill_report ( self ):
        for child in self.report_box.get_children ( ): self.report_box.remove ( child )
        report = self.backend.report ( )
        level, title, text = report["verdict"]
        banner = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 2 )
        head = Gtk.Label ( xalign = 0 )
        head.set_markup ( '<span size="large" weight="bold" foreground="{}">● {}</span>'.format ( LEVEL_COLORS[level], GLib.markup_escape_text ( title ) ) )
        body = Gtk.Label ( label = text, xalign = 0, wrap = True )
        banner.pack_start ( head, False, False, 0 )
        banner.pack_start ( body, False, False, 0 )
        self.report_box.pack_start ( banner, False, False, 0 )

        for name, levels, recommendations in report["molecules"]:
            frame = Gtk.Frame ( )
            box = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 4, margin = 8 )
            frame.add ( box )
            total = sum ( levels.values ( ) ) or 1
            head = Gtk.Label ( xalign = 0 )
            head.set_markup ( "<b>{}</b>   {}".format ( GLib.markup_escape_text ( name ), "   ".join (
                '<span foreground="{}">●</span> {} {}'.format ( LEVEL_COLORS[l], levels.get ( l, 0 ), t )
                for l, t in ( ( "green", "atoms ok" ), ( "yellow", "to review" ), ( "red", "attention" ) ) ) ) )
            box.pack_start ( head, False, False, 0 )
            bar = Gtk.LevelBar ( )
            bar.set_min_value ( 0 ); bar.set_max_value ( 1 ); bar.set_value ( levels.get ( "green", 0 ) / float ( total ) )
            bar.set_tooltip_text ( "{:.0f} % of the atoms have only database parameters".format ( 100.0 * levels.get ( "green", 0 ) / total ) )
            box.pack_start ( bar, False, False, 0 )
            for r in recommendations:
                label = Gtk.Label ( label = "• " + r, xalign = 0, wrap = True )
                box.pack_start ( label, False, False, 0 )
            self.report_box.pack_start ( frame, False, False, 0 )

        if report["issues"]:
            head = Gtk.Label ( xalign = 0 )
            head.set_markup ( "<b>Problems</b>  <small>(click to open the parameter)</small>" )
            self.report_box.pack_start ( head, False, False, 0 )
            listbox = Gtk.ListBox ( )
            listbox.set_selection_mode ( Gtk.SelectionMode.NONE )
            listbox.set_activate_on_single_click ( True )
            for level, text, term, rid in report["issues"]:
                row = Gtk.ListBoxRow ( )
                row.target = ( term, rid )
                box = Gtk.Box ( spacing = 8, margin = 4 )
                dot = Gtk.Label ( )
                dot.set_markup ( '<span foreground="{}">●</span>'.format ( LEVEL_COLORS[level] ) )
                label = Gtk.Label ( label = text, xalign = 0 )
                label.set_ellipsize ( Pango.EllipsizeMode.END )
                box.pack_start ( dot, False, False, 0 )
                box.pack_start ( label, True, True, 0 )
                box.pack_end ( Gtk.Image.new_from_icon_name ( "go-next-symbolic", Gtk.IconSize.MENU ), False, False, 0 )
                row.add ( box )
                listbox.add ( row )
            listbox.connect ( "row-activated", self.on_issue_activated )
            self.report_box.pack_start ( listbox, False, False, 0 )
        self.report_box.show_all ( )

    #-------------------------------------------------------------------------
    #  N A V I G A T I O N
    #-------------------------------------------------------------------------
    def on_sidebar_selected ( self, listbox, row ):
        if row is None or self.backend is None: return
        self.stack.set_visible_child_name ( row.name )
        is_table = row.name != "report"
        self.filter_row.set_visible ( is_table )
        self.detail.set_visible ( is_table )
        if is_table:
            self._update_filter_info ( )
            self._show_detail ( None )
            selection = self.trees[row.name].get_selection ( )
            model, tree_iter = selection.get_selected ( )
            if tree_iter is not None: self.on_row_selected ( selection, row.name )

    def on_issue_activated ( self, listbox, row ):
        term, rid = row.target
        if rid is not None:
            self.filter_buttons["all"].set_active ( True )
        else:
            self.filter_buttons["attention"].set_active ( True )
        for side in self.sidebar.get_children ( ):
            if side.name == term: self.sidebar.select_row ( side )
        if rid is None: return
        model = self.filtered[term]
        for path_row in model:
            if path_row[0] == rid:
                tree = self.trees[term]
                tree.get_selection ( ).select_iter ( path_row.iter )
                tree.scroll_to_cell ( path_row.path, None, True, 0.3, 0.0 )
                break

    def _filter_key ( self ):
        for key, button in self.filter_buttons.items ( ):
            if button.get_active ( ): return key
        return "all"

    def _row_visible ( self, model, tree_iter, data = None ):
        key = self._filter_key ( )
        n = model.get_n_columns ( )
        level, in_molecule = model[tree_iter][n - 2], model[tree_iter][n - 1]
        if key == "ligand":    return in_molecule
        if key == "attention": return level in ( "yellow", "red", "blue" )
        return True

    def on_filter_toggled ( self, button, key ):
        if not button.get_active ( ): return
        for model in self.filtered.values ( ): model.refilter ( )
        self._update_filter_info ( )

    def _update_filter_info ( self ):
        row = self.sidebar.get_selected_row ( )
        if row is None or row.name == "report" or row.name not in self.filtered:
            self.label_filter_info.set_text ( "" ); return
        shown, total = len ( self.filtered[row.name] ), len ( self.stores[row.name] )
        text = "{} of {} shown".format ( shown, total )
        if shown == 0 and self._filter_key ( ) == "attention":
            text = "No problems here -- choose All to see every parameter"
        self.label_filter_info.set_text ( text )

    #-------------------------------------------------------------------------
    #  D E T A I L   P A N E
    #-------------------------------------------------------------------------
    def on_row_selected ( self, selection, term ):
        model, tree_iter = selection.get_selected ( )
        if tree_iter is None or self.backend is None: return
        rid = model[tree_iter][0]
        if self.current and self.current[0] == term and self.current[1] == rid: return
        self.current = [ term, rid, 0 ]
        self._show_detail ( self.current )
        # . after the click is handled: centring the camera is a blocking animation
        GLib.idle_add ( lambda: ( self._project ( ), False )[1] )

    def _show_detail ( self, current ):
        for child in self.detail_fields.get_children ( ): self.detail_fields.remove ( child )
        self.button_fit.set_visible ( False )
        if current is None:
            self.detail_dot.set_text ( "" )
            self.detail_title.set_markup ( '<span foreground="#888888">Select a parameter to see where it is used and edit it.</span>' )
            self.detail_quality.set_text ( "" ); self.label_instance.set_text ( "" )
            for w in ( self.button_prev, self.button_next, self.check_show ): w.set_sensitive ( False )
            return
        term, rid, k = current
        b = self.backend
        row = b.rows ( term )[rid]
        level = b.level ( term, row )
        self.detail_dot.set_markup ( '<span foreground="{}">●</span>'.format ( LEVEL_COLORS[level] ) )
        self.detail_title.set_markup ( "<b>{}</b>".format ( GLib.markup_escape_text ( b.title ( term, row ) ) ) )
        self.detail_quality.set_text ( " · ".join ( t for t in ( b.quality_text ( term, row ), b.source_text ( term, row ) ) if t ) )
        editable = getattr ( b, "editable", True )
        for slot, label, text, unit in b.fields ( term, row ):
            self.detail_fields.pack_start ( Gtk.Label ( label = label ), False, False, 0 )
            if not editable:
                value = Gtk.Label ( xalign = 0, selectable = True )
                value.set_markup ( "<b>{}</b>".format ( GLib.markup_escape_text ( text ) ) )
                self.detail_fields.pack_start ( value, False, False, 0 )
                continue
            entry = Gtk.Entry ( text = text )
            entry.set_width_chars ( 26 if term == "dihedral" else 10 )
            entry.set_tooltip_text ( "Type a new value and press Enter; Apply writes the changes to the system" )
            entry.connect ( "activate", self.on_field_edited, term, rid, slot )
            entry.connect ( "focus-out-event", lambda w, e, t = term, r = rid, s = slot: ( self.on_field_edited ( w, t, r, s ), False )[1] )
            self.detail_fields.pack_start ( entry, False, False, 0 )
            unit_label = Gtk.Label ( label = unit )
            unit_label.get_style_context ( ).add_class ( "dim-label" )
            self.detail_fields.pack_start ( unit_label, False, False, 0 )
        self.detail_fields.show_all ( )
        n = len ( row.get ( "instances", [ ] ) )
        self.button_prev.set_sensitive ( n > 1 ); self.button_next.set_sensitive ( n > 1 )
        self.check_show.set_sensitive ( n > 0 )
        self.button_fit.set_visible ( b.can_fit_torsion ( term ) )

    def on_field_edited ( self, entry, term, rid, slot ):
        if self.backend is None or not self.current or self.current[:2] != [ term, rid ]: return
        row = self.backend.rows ( term )[rid]
        try:
            changed = self.backend.set_field ( term, row, slot, entry.get_text ( ).strip ( ) )
        except ( ValueError, IndexError ):
            self.label_status.set_markup ( '<span foreground="{}">Invalid value "{}".</span>'.format (
                LEVEL_COLORS["red"], GLib.markup_escape_text ( entry.get_text ( ) ) ) )
            return
        if not changed: return
        self._refresh_row ( term, rid )
        self._update_badges ( )
        self._update_footer ( )
        self._project ( )

    def _refresh_row ( self, term, rid ):
        b = self.backend
        row = b.rows ( term )[rid]
        level = b.level ( term, row )
        for store_row in self.stores[term]:
            if store_row[0] == rid:
                values = [ rid, '<span foreground="{}">●</span>'.format ( LEVEL_COLORS[level] ) ] + \
                         [ func ( row ) for title, func in b.columns ( term ) ] + [ b.quality_text ( term, row ), level, b.in_molecule ( term, row ) ]
                for k, v in enumerate ( values ): store_row[k] = v
        self.detail_dot.set_markup ( '<span foreground="{}">●</span>'.format ( LEVEL_COLORS[level] ) )
        self.detail_quality.set_text ( "{} · {}".format ( b.quality_text ( term, row ), b.source_text ( term, row ) ) )

    def _update_footer ( self ):
        b = self.backend
        n = b.n_pending ( ) if b is not None else 0
        self.button_apply.set_sensitive ( n > 0 )
        self.button_revert.set_sensitive ( n > 0 )
        if b is None:
            self.label_status.set_text ( "" ); return
        total = b.total_charge ( ) + 0.0
        if abs ( total ) < 5.0e-5: total = 0.0
        integral = abs ( total - round ( total ) ) <= 1.0e-4
        charge = '<span foreground="{}">total charge {:+.4f}</span>'.format ( LEVEL_COLORS["green"] if integral else LEVEL_COLORS["red"], total )
        if n:
            text = "{} change(s) not applied · {}".format ( n, charge )
            if not integral: text += " -- must be an integer before Apply"
        else:
            text = charge
        self.label_status.set_markup ( text )

    #-------------------------------------------------------------------------
    #  3 D   V I E W
    #-------------------------------------------------------------------------
    def _vobject ( self ):
        e_id = getattr ( self.backend.system, "e_id", None )
        candidates = [ v for v in self.vm_session.vm_objects_dic.values ( ) if getattr ( v, "e_id", None ) == e_id ]
        if not candidates: return None
        active = [ v for v in candidates if getattr ( v, "active", False ) ]
        return ( active or candidates )[-1]

    def _step ( self, delta ):
        if not self.current: return
        term, rid, k = self.current
        n = len ( self.backend.rows ( term )[rid].get ( "instances", [ ] ) ) or 1
        self.current[2] = ( k + delta ) % n
        self._project ( )

    def _project ( self ):
        try:
            self._project_inner ( )
        except Exception:
            traceback.print_exc ( )

    def _project_inner ( self ):
        import numpy as np
        if not self.current or self.backend is None: return
        term, rid, k = self.current
        row = self.backend.rows ( term )[rid]
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
        text = "{}/{}  {}{}".format ( k + 1, len ( instances ), residue, "-".join ( names ) )
        if atoms and all ( a is not None for a in atoms ):
            points = [ a.coords ( ) for a in atoms ]
            measured = self.backend.measure_text ( term, row, points )
            if measured: text += "   ·   " + measured
            if self.check_show.get_active ( ) and term != "lj":
                self._pick ( atoms, vobject, points )
        self.label_instance.set_text ( text )

    def _pick ( self, atoms, vobject, points ):
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
        centre = np.mean ( np.array ( points, dtype = np.float32 ), axis = 0 ).astype ( np.float32 )
        try:
            session.vm_glcore.center_on_coordinates ( vobject, centre )
        except Exception:
            traceback.print_exc ( )
        session.vm_glcore.queue_draw ( )

    #-------------------------------------------------------------------------
    #  A C T I O N S
    #-------------------------------------------------------------------------
    def on_revert ( self, button ):
        self.on_system_changed ( )
        self.label_status.set_text ( "Changes discarded." )

    def on_apply ( self, button ):
        b = self.backend
        if b is None or not b.n_pending ( ): return
        system = b.system
        try:
            before = system.Energy ( log = None )
        except Exception:
            before = None
        try:
            summary = b.apply ( )
        except Exception as error:
            traceback.print_exc ( )
            self.main.simple_dialog.error ( msg = str ( error ) )
            return
        try:
            after = system.Energy ( log = None )
        except Exception:
            after = None
        page = self.sidebar.get_selected_row ( )
        page = page.name if page is not None else "report"
        self.on_system_changed ( )
        for row in self.sidebar.get_children ( ):
            if row.name == page: self.sidebar.select_row ( row )
        energy = "  Energy {:.2f} -> {:.2f} kJ/mol.".format ( before, after ) if before is not None and after is not None else ""
        self.label_status.set_text ( "Applied: {}.{}".format ( summary, energy ) )
        try:
            self.main.bottom_notebook.status_teeview_add_new_item (
                message = "{} parameters edited: {}.".format ( b.force_field, summary ), system = system )
            self.main.refresh_main_statusbar ( )
        except Exception:
            pass

    def on_fit_torsion ( self, button ):
        from gui.windows.builder.torsion_fit_window import TorsionFitWindow
        import numpy as np
        if not self.current or self.backend is None or not self.backend.can_fit_torsion ( self.current[0] ): return
        if self.backend.n_pending ( ):
            self.label_status.set_text ( "Apply or revert the pending changes before fitting a torsion." ); return
        term, rid, k = self.current
        row = self.backend.rows ( term )[rid]
        coordinates = None
        vobject = self._vobject ( )
        if vobject is not None and len ( vobject.atoms ) == len ( self.backend.system.atoms ):
            frame = min ( getattr ( self.vm_session, "frame", 0 ) or 0, len ( vobject.frames ) - 1 )
            coordinates = np.array ( vobject.frames[frame], dtype = float )
        try:
            target = self.backend.torsion_target ( row, row["instances"][k], coordinates )
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
            self.main.bottom_notebook.status_teeview_add_new_item ( message = "OPLS " + message, system = getattr ( self.backend, "system", None ) )
            self.main.refresh_main_statusbar ( )
        except Exception:
            pass

    def _open_path ( self, path ):
        if not path:
            self.label_status.set_text ( "Not available for this system." ); return
        try:
            subprocess.Popen ( [ "xdg-open", path ] )
        except Exception as error:
            self.label_status.set_text ( "Could not open {}: {}".format ( path, error ) )

    def on_open_report ( self, item ):
        self._open_path ( self.backend.report_file ( ) if self.backend else None )

    def on_open_folder ( self, item ):
        self._open_path ( self.backend.folder ( ) if self.backend else None )
