#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  periodic_table_dialog.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Full periodic table element picker for the Builder (user's own request:
#  "um botao onde podemos selecionar qualquer elemento quimico" -- the
#  sidebar's 4 quick-pick radios, C/N/O/H, stay as they are, this is the
#  "anything else" escape hatch, matching the Builder roadmap's own
#  Section 7 note: "keeping C/N/O/H as one-click quick buttons").
#
#  [EN] "use os parametros do UFF como referencia" (user's own explicit
#  request): confirmed by direct research (not guessed) that this
#  project's own vismol.utils.elements.PeriodicTable ALREADY carries
#  UFF-derived per-element data for the WHOLE periodic table -- columns
#  7/8 of every entry ("r_UFF"/"En_UFF") are copied from Rappe et al.'s
#  original 1992 UFF paper, the SAME source pDynamo3's own DYFF force
#  field's atomTypes.yaml cites verbatim ("Atom types copied from UFF and
#  put into DYFF-1.0 format" -- parameters/forceFields/dyff/dyff-1.0/
#  atomTypes.yaml:1). Atom.cov_rad/vdw_rad/electronegativity (vismol/
#  model/atom.py's _init_cov_rad()/_init_vdw_rad()/_init_electronegativity())
#  already read straight from this SAME table for every atom this Builder
#  ever creates -- so bond-length rendering and DYFF assignment already
#  work correctly for ANY element picked here, no new per-element data
#  needed for that part. The only place this project's OWN element
#  coverage was ever narrower than UFF's is atom_ops.STANDARD_VALENCE
#  (used for automatic hydrogen add/remove) -- extended separately, in
#  atom_ops.py itself, using the same UFF-derived "typical" valence for
#  the additional main-group elements added there.
#
#  Capped at Z=103 (Lr) -- matches pDynamo3's own DYFF/UFF parameter
#  set's documented coverage (dyff-1.0/atomTypes.yaml goes up to Z=103)
#  exactly, so every element this picker offers is guaranteed to also be
#  typeable by DYFF later (see project_dyff_force_field_assignment) --
#  the handful of post-Lr synthetic superheavy elements have no
#  established chemistry to speak of anyway.
# ============================================================================
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk, Gdk


# [EN] Standard periodic table SHAPE (period/row -> ordered (symbol, group
# column 1-18) pairs) -- this is fixed, universally standardised chemistry
# data (not something to "look up" per app), hand-encoded once here purely
# to lay the picker out recognisably, same visual convention as any real
# periodic table wall chart. The lanthanide/actinide rows are drawn
# separately below the main 7x18 grid, same convention.
_MAIN_ROWS = [
    [ ( "H", 1 ), ( "He", 18 ) ],
    [ ( "Li", 1 ), ( "Be", 2 ), ( "B", 13 ), ( "C", 14 ), ( "N", 15 ), ( "O", 16 ), ( "F", 17 ), ( "Ne", 18 ) ],
    [ ( "Na", 1 ), ( "Mg", 2 ), ( "Al", 13 ), ( "Si", 14 ), ( "P", 15 ), ( "S", 16 ), ( "Cl", 17 ), ( "Ar", 18 ) ],
    [ ( "K", 1 ), ( "Ca", 2 ), ( "Sc", 3 ), ( "Ti", 4 ), ( "V", 5 ), ( "Cr", 6 ), ( "Mn", 7 ), ( "Fe", 8 ),
      ( "Co", 9 ), ( "Ni", 10 ), ( "Cu", 11 ), ( "Zn", 12 ), ( "Ga", 13 ), ( "Ge", 14 ), ( "As", 15 ),
      ( "Se", 16 ), ( "Br", 17 ), ( "Kr", 18 ) ],
    [ ( "Rb", 1 ), ( "Sr", 2 ), ( "Y", 3 ), ( "Zr", 4 ), ( "Nb", 5 ), ( "Mo", 6 ), ( "Tc", 7 ), ( "Ru", 8 ),
      ( "Rh", 9 ), ( "Pd", 10 ), ( "Ag", 11 ), ( "Cd", 12 ), ( "In", 13 ), ( "Sn", 14 ), ( "Sb", 15 ),
      ( "Te", 16 ), ( "I", 17 ), ( "Xe", 18 ) ],
    [ ( "Cs", 1 ), ( "Ba", 2 ), ( "La", 3 ), ( "Hf", 4 ), ( "Ta", 5 ), ( "W", 6 ), ( "Re", 7 ), ( "Os", 8 ),
      ( "Ir", 9 ), ( "Pt", 10 ), ( "Au", 11 ), ( "Hg", 12 ), ( "Tl", 13 ), ( "Pb", 14 ), ( "Bi", 15 ),
      ( "Po", 16 ), ( "At", 17 ), ( "Rn", 18 ) ],
    [ ( "Fr", 1 ), ( "Ra", 2 ), ( "Ac", 3 ) ],   # capped at Z=103 -- see module docstring
]
_LANTHANIDES = [ "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu" ]
_ACTINIDES   = [ "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr" ]
_F_BLOCK_FIRST_COLUMN = 4   # matches the real periodic table's own footnote-row convention

# [EN] Quick-pick elements (already have their OWN dedicated sidebar
# radio buttons -- builder_sidebar.glade's element_box) get a visually
# distinct border here, so it's obvious picking them here is redundant
# with, not different from, those buttons.
_QUICK_PICK_SYMBOLS = { "C", "N", "O", "H" }


def _element_tile_color ( vm_session, symbol ):
    """ Same per-element RGB this Builder already uses to COLOUR that
    element's own atoms in the 3D view (vm_session.periodic_table --
    see this module's own docstring for where that table's data comes
    from) -- so a tile's colour here is a preview of how it'll actually
    look once placed, not an arbitrary picker-only colour scheme. """
    entry = vm_session.periodic_table.elements_by_symbol.get ( symbol )
    if entry is None:
        return ( 0.7, 0.7, 0.7 )
    return tuple ( entry[2] )


def choose_element ( main, parent_window = None ):
    """ Shows the periodic table picker modally and returns the chosen
    element symbol (e.g. "Fe"), or None if the user cancelled. Doesn't
    itself touch vm_session.builder_atom_symbol -- see builder_sidebar.
    py's own on_periodic_table_button_clicked() for that, so this
    function stays reusable anywhere else an element needs picking. """
    vm_session = main.vm_session

    dialog = Gtk.Dialog ( title = "Periodic Table", transient_for = parent_window or main.window, modal = True )
    dialog.add_buttons ( Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL )
    dialog.set_default_size ( 560, 320 )

    content = dialog.get_content_area ( )
    content.set_spacing ( 6 )
    content.set_border_width ( 8 )

    # --- search row: type a symbol or (part of) a name, Enter jumps
    # straight to it if there's exactly one match, same convenience a
    # much bigger grid (118 tiny buttons) badly needs. ---
    search_box = Gtk.Box ( orientation = Gtk.Orientation.HORIZONTAL, spacing = 6 )
    search_entry = Gtk.Entry ( )
    search_entry.set_placeholder_text ( "Search by symbol or name (e.g. \"Fe\" or \"iron\") -- Enter to pick" )
    search_box.pack_start ( search_entry, True, True, 0 )
    content.pack_start ( search_box, False, False, 0 )

    grid = Gtk.Grid ( )
    grid.set_row_spacing ( 2 )
    grid.set_column_spacing ( 2 )
    grid.set_row_homogeneous ( True )
    grid.set_column_homogeneous ( True )

    scroller = Gtk.ScrolledWindow ( )
    scroller.set_policy ( Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC )
    scroller.add ( grid )
    content.pack_start ( scroller, True, True, 0 )

    result = { "symbol": None }

    def _pick ( symbol ):
        result["symbol"] = symbol
        dialog.response ( Gtk.ResponseType.OK )

    def _add_tile ( symbol, row, col ):
        entry = vm_session.periodic_table.elements_by_symbol.get ( symbol )
        name = entry[1] if entry is not None else symbol
        number = entry[0] if entry is not None else "?"

        button = Gtk.Button ( label = symbol )
        button.set_tooltip_text ( "{} ({}) -- Z={}".format ( name, symbol, number ) )
        r, g, b = _element_tile_color ( vm_session, symbol )
        css = Gtk.CssProvider ( )
        border = "2px solid #333333" if symbol in _QUICK_PICK_SYMBOLS else "1px solid #999999"
        css.load_from_data ( "button {{ background-color: rgba({},{},{},0.55); border: {}; padding: 2px; }}".format (
                int ( r * 255 ), int ( g * 255 ), int ( b * 255 ), border ).encode ( "utf-8" ) )
        button.get_style_context ( ).add_provider ( css, Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION )
        button.connect ( "clicked", lambda w, s = symbol: _pick ( s ) )
        grid.attach ( button, col - 1, row - 1, 1, 1 )

    for row_index, row_entries in enumerate ( _MAIN_ROWS, start = 1 ):
        for ( symbol, col ) in row_entries:
            _add_tile ( symbol, row_index, col )

    lanthanide_row = len ( _MAIN_ROWS ) + 2   # +2, not +1: leaves an empty spacer row, matching every real periodic table chart
    actinide_row   = lanthanide_row + 1
    for offset, symbol in enumerate ( _LANTHANIDES ):
        _add_tile ( symbol, lanthanide_row, _F_BLOCK_FIRST_COLUMN + offset )
    for offset, symbol in enumerate ( _ACTINIDES ):
        _add_tile ( symbol, actinide_row, _F_BLOCK_FIRST_COLUMN + offset )

    def _on_search_activate ( entry ):
        query = entry.get_text ( ).strip ( ).lower ( )
        if not query:
            return
        by_symbol = vm_session.periodic_table.elements_by_symbol
        exact_symbol = [ s for s in by_symbol.keys ( ) if s.lower ( ) == query ]
        if exact_symbol:
            _pick ( exact_symbol[0] )
            return
        name_matches = [ s for s, e in by_symbol.items ( ) if query in e[1].lower ( ) ]
        if len ( name_matches ) == 1:
            _pick ( name_matches[0] )
        elif not name_matches:
            entry.get_style_context ( ).add_class ( "error" )
        else:
            entry.get_style_context ( ).add_class ( "error" )   # ambiguous -- multiple names contain the query
    search_entry.connect ( "activate", _on_search_activate )

    dialog.show_all ( )
    response = dialog.run ( )
    dialog.destroy ( )

    if response == Gtk.ResponseType.OK:
        return result["symbol"]
    return None
