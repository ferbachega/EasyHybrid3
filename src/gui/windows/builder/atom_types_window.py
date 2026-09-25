#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  atom_types_window.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  "Atom Types" window -- user's own explicit request: for a SELECTION
#  of atoms (the app's own existing "viewing selection" mechanism --
#  shift+click / shift+drag, works alongside any active Builder tool
#  with no conflict) inside the Builder, show which DYFF atom type each
#  one currently perceives as, and let the user CORRECT it by hand,
#  directly in an editable treeview cell, if the automatic perception
#  got it wrong (the guanidinium nitrogen mistyping this project already
#  found+fixed with a hardcoded structural detector is exactly the class
#  of mistake this window lets the user fix themselves, for ANY atom,
#  without needing a new one-off detector written for every future case).
#
#  Mirrors the exact open_window()/close_window() self-contained Glade+
#  Python window pattern already used by builder_sidebar.py/fragment_
#  library_window.py.
# ============================================================================
import os
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk

from gui.windows.builder.atom_types import (
    compute_atom_types, set_manual_atom_type_override, valid_type_labels_for_symbol,
)


class AtomTypesWindow ( ):
    """ Atom Types inspector/editor. """

    def __init__ ( self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.visible    = False
        self.target_object = None
        self.atom_ids      = [ ]

    def open_window ( self, target_object, atom_ids ):
        """ (Re)opens the window showing `atom_ids` (a list/set of
        atom_id ints, ALREADY filtered by the caller to atoms belonging
        to `target_object` -- see builder_sidebar.py's own
        on_atom_types_button_clicked(), which does that filtering
        against the app's existing viewing-selection mechanism). Re-
        populates from scratch on every call (even if already visible,
        same "re-scan rather than trust stale state" philosophy
        fragment_library_window.py's own open_window() already uses),
        so re-selecting a different set of atoms and clicking the
        sidebar button again always shows the CURRENT selection, not
        whatever was open before.

        [EN] 2026-09-25: stores the ATOM OBJECTS themselves (self.atoms),
        not just their atom_id ints -- protonate_atom()/deprotonate_
        atom() (see on_protonate_button_clicked()/_deprotonate_button_
        clicked()) can REMOVE a hydrogen, which renumbers every atom_id
        ABOVE it down by one (atom_ops.remove_atom()'s own docstring) --
        but keeps the SAME Python object for every surviving atom, only
        updating that object's OWN `.atom_id` attribute in place. Re-
        deriving atom_ids from these objects' current `.atom_id` on every
        _refresh_rows() call (instead of trusting a list of ints captured
        once, here) keeps every row pointing at the atom the user
        actually selected, even after a shift earlier in the same
        multi-row batch. """
        self.target_object = target_object
        self.atoms = [ target_object.atoms[atom_id] for atom_id in sorted ( atom_ids )
                        if atom_id in target_object.atoms ]
        self.atom_ids = [ ]   # populated fresh by _refresh_rows() below

        if not self.visible:
            self.builder = Gtk.Builder ( )
            self.builder.add_from_file ( os.path.join ( self.main.home, 'src/gui/windows/builder/atom_types_window.glade' ) )
            self.builder.connect_signals ( self )

            self.window        = self.builder.get_object ( 'window' )
            self.liststore      = self.builder.get_object ( 'atom_types_liststore' )
            self.treeview        = self.builder.get_object ( 'atom_types_treeview' )
            self.treeview_selection = self.builder.get_object ( 'atom_types_treeview_selection' )
            self.header_label   = self.builder.get_object ( 'header_label' )
            self.status_label   = self.builder.get_object ( 'status_label' )

            self.window.show_all ( )
            self.visible = True
        else:
            self.window.present ( )

        self._refresh_rows ( )

    def _refresh_rows ( self ):
        """ Re-runs compute_atom_types() (fresh DYFF perception straight
        off the target object's CURRENT structure -- see atom_types.py's
        own docstring for why this doesn't need a prior real DYFF
        assignment) and repopulates every row from scratch.

        [EN] 2026-09-25: self.atom_ids is derived FRESH here from self.
        atoms's own current `.atom_id` (see open_window()'s own comment)
        -- an atom whose id now maps to a DIFFERENT object (this one was
        actually deleted by something else, e.g. the Delete tool, and its
        old id slot got reused by a later atom shifting down) is dropped,
        same as the pre-existing "deleted by some other action" handling
        further down already did for a stale int id. """
        if self.target_object is None:
            return

        self.atom_ids = [ atom.atom_id for atom in self.atoms
                           if self.target_object.atoms.get ( atom.atom_id ) is atom ]

        self.header_label.set_text ( "{} -- {} atom(s) selected".format (
                self.target_object.name, len ( self.atom_ids ) ) )

        types_by_id = compute_atom_types ( self.target_object, atom_ids = self.atom_ids )
        overrides   = getattr ( self.target_object, "manual_atom_type_overrides", None ) or { }

        # [EN] User's own follow-up request -- "seria interessante colocar
        # o label dos indexes no vobject, mas [para] serem facilmente
        # identificados": while this window is open, the SAME atom_id
        # shown in column 0 here is ALSO drawn as a numbered ring+label
        # right on top of that atom in the 3D view (vm_session.
        # builder_atom_types_labeled_atoms, read by vismol_glcore.py's
        # render() -- reuses click_mode.draw_hover_highlight()/
        # draw_ring_label(), the EXACT same overlay mechanism the Bond
        # Order/Rotate Dihedral tools already use for their own armed-
        # atom previews), so a row in this table and the atom it refers
        # to on screen are trivially easy to match up, instead of having
        # to guess or click each one to check.
        labeled_atoms = [ ]

        self.liststore.clear ( )
        for atom_id in self.atom_ids:
            atom = self.target_object.atoms.get ( atom_id )
            if atom is None:
                continue   # deleted by some other action since the selection was made
            labeled_atoms.append ( atom )
            perceived, effective = types_by_id.get ( atom_id, ( None, None ) )
            is_override = atom_id in overrides
            valid_labels = valid_type_labels_for_symbol ( self.vm_session, atom.symbol )
            tooltip = "Valid {} types: {}".format ( atom.symbol, ", ".join ( valid_labels ) ) if valid_labels \
                      else "No DYFF type is defined for element {}.".format ( atom.symbol )
            self.liststore.append ( [
                    atom_id,
                    atom.symbol,
                    atom.name,
                    effective if effective is not None else "?",
                    "manual" if is_override else "auto",
                    tooltip,
                    int ( getattr ( atom, "formal_charge", 0 ) ),
            ] )

        self.vm_session.builder_atom_types_labeled_atoms = labeled_atoms
        if getattr ( self.vm_session, "vm_glcore", None ) is not None:
            self.vm_session.vm_glcore.queue_draw ( )

    # ------------------------------------------------------------------
    #  Signal handlers (referenced by name in atom_types_window.glade)
    # ------------------------------------------------------------------

    def on_treeview_selection_changed ( self, selection ):
        """ [EN] User's own explicit request -- "ao clicar sobre um
        atomo, ele deve automaticamente centralizar no atomo (na
        GLArea), isso vai ajudar a identificar o tipo atomico e o
        ambiente": clicking a row moves the 3D camera to centre on that
        SAME atom, reusing vm_glcore.center_on_atom() -- the EXACT same
        primitive the pre-existing "Go To Atom" window already uses for
        its own atom treeview (gui/windows/setup/windows_and_dialogs/
        selection_windows/go_to_atom.py), not a new camera-centring
        mechanism.

        Centres on the LAST row in the current selection (get_selected_
        rows(), not get_selected() -- this treeview's own selection mode
        is "multiple", for "Clear Override(s)"'s own multi-row support,
        and get_selected() only ever works reliably in SINGLE selection
        mode) -- for a plain click that's simply the row just clicked;
        for a shift/ctrl-extended multi-select it's whichever row the
        selection was just extended to, which is the one the user's
        cursor/attention is actually on. """
        model, paths = selection.get_selected_rows ( )
        if not paths or self.target_object is None:
            return
        atom_id = model[paths[-1]][0]
        atom = self.target_object.atoms.get ( atom_id )
        if atom is None:
            return
        vm_glcore = getattr ( self.vm_session, "vm_glcore", None )
        if vm_glcore is not None:
            vm_glcore.center_on_atom ( atom )

    def on_type_cell_edited ( self, renderer, path, new_text ):
        """ Validates new_text against this row's OWN element's valid
        DYFF labels (valid_type_labels_for_symbol()) before accepting
        it -- an invalid label is rejected with a clear message in the
        status label and the cell is left UNCHANGED (neither the
        liststore nor manual_atom_type_overrides are touched), rather
        than silently storing something DYFF's own typer would later
        choke on. Clearing the cell (empty text) removes any existing
        override for that atom instead, reverting to the auto-perceived
        type. """
        new_text = new_text.strip ( )
        atom_id  = self.liststore[path][0]
        atom     = self.target_object.atoms.get ( atom_id ) if self.target_object else None
        if atom is None:
            return

        if not new_text:
            set_manual_atom_type_override ( self.target_object, atom_id, None )
            self.status_label.set_text ( "Atom #{}: override cleared -- back to the auto-perceived type.".format ( atom_id ) )
            self._refresh_rows ( )
            return

        valid_labels = valid_type_labels_for_symbol ( self.vm_session, atom.symbol )
        if valid_labels and new_text not in valid_labels:
            self.status_label.set_text (
                    "\"{}\" is not a valid DYFF type for {} -- valid options: {}.".format (
                            new_text, atom.symbol, ", ".join ( valid_labels ) ) )
            return

        set_manual_atom_type_override ( self.target_object, atom_id, new_text )
        self.status_label.set_text ( "Atom #{}: type manually set to \"{}\".".format ( atom_id, new_text ) )
        self._refresh_rows ( )

    def on_charge_cell_edited ( self, renderer, path, new_text ):
        """ [EN] 2026-09-25, user's own explicit request: "adicionar a
        possibilidade de alterarmos estado de protonacao". A RAW override
        -- sets atom.formal_charge directly, WITHOUT touching hydrogens
        (unlike the Protonate/Deprotonate buttons below, which change
        both together) -- same "manual correction, no automatic side
        effects" philosophy as on_type_cell_edited() just above. Useful
        for fixing a charge on an atom whose element isn't in atom_ops.
        STANDARD_VALENCE at all (Protonate/Deprotonate can't touch those),
        or for setting a charge with no accompanying H change. Rejects
        non-integer text, leaving the cell/atom unchanged. """
        new_text = new_text.strip ( )
        atom_id  = self.liststore[path][0]
        atom     = self.target_object.atoms.get ( atom_id ) if self.target_object else None
        if atom is None:
            return
        try:
            new_charge = int ( new_text )
        except ValueError:
            self.status_label.set_text ( "\"{}\" is not a valid integer charge.".format ( new_text ) )
            return

        atom.formal_charge = new_charge
        from gui.windows.builder.empty_object import sync_pdynamo_system
        sync_pdynamo_system ( self.target_object )
        self.status_label.set_text ( "Atom #{}: formal charge set to {:+d}.".format ( atom_id, new_charge ) )
        self._refresh_rows ( )

    def on_refresh_button_clicked ( self, button ):
        self._refresh_rows ( )
        self.status_label.set_text ( "Refreshed." )

    def on_clear_override_button_clicked ( self, button ):
        """ Clears the override for every currently SELECTED row (treeview
        selection, not the Builder's own viewing selection) -- lets the
        user fix several rows at once instead of one cell at a time. """
        model, paths = self.treeview_selection.get_selected_rows ( )
        cleared = 0
        for path in paths:
            atom_id = model[path][0]
            set_manual_atom_type_override ( self.target_object, atom_id, None )
            cleared += 1
        if cleared:
            self.status_label.set_text ( "Cleared {} override(s).".format ( cleared ) )
            self._refresh_rows ( )
        else:
            self.status_label.set_text ( "Select one or more rows first." )

    def _protonation_atoms_from_selection ( self ):
        """ [EN] Shared by on_protonate_button_clicked()/on_deprotonate_
        button_clicked() below: resolves the treeview's CURRENTLY
        selected row(s) to actual Atom OBJECTS (not bare atom_id ints).
        This matters because atom_ops.deprotonate_atom() can REMOVE a
        hydrogen with a LOWER atom_id than another still-pending selected
        atom -- which shifts every atom_id above it down by one (see
        atom_ops.remove_atom()'s own docstring) -- but the SAME Atom
        object stays the same Python instance throughout (only its own
        `.atom_id` attribute gets updated in place), so re-reading
        `atom.atom_id` fresh right before each individual protonate/
        deprotonate call (done by the callers below, not here) always
        targets the correct, still-live atom regardless of how many
        earlier removals in this same batch shifted ids around it. """
        model, paths = self.treeview_selection.get_selected_rows ( )
        atoms = [ ]
        for path in paths:
            atom_id = model[path][0]
            atom = self.target_object.atoms.get ( atom_id ) if self.target_object else None
            if atom is not None:
                atoms.append ( atom )
        return atoms

    def on_protonate_button_clicked ( self, button ):
        """ [EN] 2026-09-25, user's own explicit request: "adicionar a
        possibilidade de alterarmos estado de protonacao, especialmente e
        aminas, que em agua deveriam aceitar um nitrogenio tetravalente e
        com carga total +1" -- calls atom_ops.protonate_atom() for every
        currently SELECTED row (treeview selection, same multi-row
        pattern as Clear Override(s) above). See _protonation_atoms_from_
        selection()'s own docstring for why each atom is re-resolved to
        its CURRENT atom_id (via the captured Atom object's own `.
        atom_id`) right before each call, instead of using the atom_id
        values captured before the loop started. """
        atoms = self._protonation_atoms_from_selection ( )
        if not atoms:
            self.status_label.set_text ( "Select one or more rows first." )
            return

        from gui.windows.builder.atom_ops import protonate_atom
        messages = [ ]
        for atom in atoms:
            ok, message = protonate_atom ( self.target_object, atom.atom_id )
            messages.append ( message )
        self.status_label.set_text ( "  |  ".join ( messages ) )
        self._refresh_rows ( )

    def on_deprotonate_button_clicked ( self, button ):
        """ Symmetric counterpart to on_protonate_button_clicked() above
        -- see that method's and atom_ops.deprotonate_atom()'s own
        docstrings. """
        atoms = self._protonation_atoms_from_selection ( )
        if not atoms:
            self.status_label.set_text ( "Select one or more rows first." )
            return

        from gui.windows.builder.atom_ops import deprotonate_atom
        messages = [ ]
        for atom in atoms:
            ok, message = deprotonate_atom ( self.target_object, atom.atom_id )
            messages.append ( message )
        self.status_label.set_text ( "  |  ".join ( messages ) )
        self._refresh_rows ( )

    def on_close_button_clicked ( self, *args ):
        """ Accepts *args so it can be connected directly to BOTH the
        "Close" button's `clicked` signal (button) -> and the window's
        own `delete-event` (widget, event) -- same convention every
        other self-contained Builder window in this project already
        uses. Deliberately does NOT clear manual_atom_type_overrides --
        closing this window is just putting it away, not discarding the
        user's own corrections. DOES clear builder_atom_types_labeled_
        atoms (the 3D-view index-label overlay, see _refresh_rows()'s
        own comment) -- those labels only make sense while this window
        is actually showing the rows they correspond to. """
        self.vm_session.builder_atom_types_labeled_atoms = [ ]
        if getattr ( self.vm_session, "vm_glcore", None ) is not None:
            self.vm_session.vm_glcore.queue_draw ( )
        if getattr ( self, "window", None ) is not None:
            self.window.destroy ( )
        self.visible = False
        return True
