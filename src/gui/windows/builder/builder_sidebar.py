#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  builder_sidebar.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Builder editing sidebar -- small utility window with the controls for
#  everything built into gui/windows/builder/{atom_ops,click_mode}.py so
#  far: turning "placemode" (vm_session.builder_atom_mode) on/off,
#  picking the active tool (Add / Delete -- vm_session.builder_tool),
#  picking the element new atoms are placed as (vm_session.
#  builder_atom_symbol), and Undo (atom_ops.undo()).
#
#  [EN] DESIGN CHOICE: opening this window IS what puts EasyHybrid into
#  Builder editing mode (matching the request that this panel's
#  PRESENCE tracks whether placemode is active, rather than being an
#  always-visible panel with its own separate open/closed state) -- if
#  there's no vm_session.builder_target_object yet, one is created
#  automatically (empty_object.create_empty_vismol_object(), the exact
#  same helper the terminal's `new` command already uses). Closing the
#  window turns editing back off. The "Editing: ON/OFF" toggle INSIDE
#  the sidebar is a separate, lighter-weight PAUSE: it stops atom-
#  placement/deletion clicks from doing anything without closing the
#  window or losing your current tool/element selection -- handy for
#  just rotating/inspecting the molecule for a moment mid-edit.
#
#  Mirrors the exact open_window()/close_window() pattern already used
#  by setup/easyhybrid_terminal.py's TerminalWindow (same Gtk.Builder +
#  add_from_file + connect_signals(self) flow), for consistency with the
#  rest of this codebase's auxiliary windows.
# ============================================================================
import os
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk

from gui.windows.builder.empty_object import create_empty_vismol_object
from gui.windows.builder.empty_object import discard_builder_object_if_unused
from gui.windows.builder.empty_object import finish_editing_existing_system
from gui.windows.builder.empty_object import hide_other_vobjects_for_builder
from gui.windows.builder.empty_object import restore_hidden_vobjects_for_builder
from gui.windows.builder.atom_ops     import undo as atom_ops_undo
from gui.windows.builder.atom_ops     import clean_up_structure
from gui.windows.builder.atom_ops     import optimize_geometry_dyff
from gui.windows.builder.fragment_library import load_fragment, FragmentError

# [EN] 2026-09-24/25: quick-pick fragment row (sidebar's own
# fragment_grid) -- maps each radio's widget attribute name to the
# ".mol2" file it loads (fragments/functional_groups/, same starter set
# the full Fragment Library browses). Only the 7 fragments confirmed
# loadable under fragment_library.py's own non-hypervalent
# _STANDARD_VALENCE table (see that module's docstring) are mapped here
# -- carboxylate/-SO2/-PO3/-NO2/-CON (this round's own "-CON" button,
# frag_con_radio1, was added to the glade with no signal handler wired
# yet either) need either a second attachment point or a charged/
# hypervalent atom this Builder doesn't support yet. 2026-09-25, user's
# own explicit request: "podemos linkar os fragmentos aos botoes em si
# mais tarde" -- fragment_grid now has MORE buttons than real mappings
# (see builder_sidebar.glade -- frag_coo_radio/frag_con_radio1/
# frag_coo_radio1/frag_coo_radio2), left deliberately unmapped/inert for
# now rather than guessed at; on_fragment_quick_changed() below already
# no-ops safely on any button with no entry here.
_QUICK_FRAGMENT_FILES = {
    "frag_cooh_radio": "carboxyl.mol2",
    "frag_nco_radio":  "isocyanate.mol2",
    "frag_ome_radio":  "methoxy.mol2",
    "frag_hex_radio":  "hexyl.mol2",
    "frag_pent_radio": "pentyl.mol2",
    "frag_benz_radio": "phenyl.mol2",
    "frag_furo_radio": "furyl.mol2",
}


class BuilderSidebarWindow ( ):
    """ Builder editing sidebar. """

    def __init__ ( self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.visible    = False

    def open_window ( self ):
        """ Opens the sidebar -- and, per this module's own design note
        above, turns Builder editing mode ON, creating a fresh empty
        target object first if none exists yet. """
        if self.visible:
            self.window.present ( )
            return

        self.builder = Gtk.Builder ( )
        self.builder.add_from_file ( os.path.join ( self.main.home, 'src/gui/windows/builder/builder_sidebar.glade' ) )
        self.builder.connect_signals ( self )

        self.window = self.builder.get_object ( 'window' )
        self.window.set_title ( 'Builder' )
        self.window.set_keep_above ( True )

        self.placemode_toggle  = self.builder.get_object ( 'placemode_toggle' )
        self.tool_add_radio    = self.builder.get_object ( 'tool_add_radio' )
        self.tool_delete_radio = self.builder.get_object ( 'tool_delete_radio' )
        self.tool_move_radio   = self.builder.get_object ( 'tool_move_radio' )
        self.tool_bond_order_radio = self.builder.get_object ( 'tool_bond_order_radio' )
        # [EN] 2026-09-25: no more dedicated "Fragment" TOOL radio -- see
        # this module's own _QUICK_FRAGMENT_FILES comment above.
        self.tool_add_structure_radio   = self.builder.get_object ( 'tool_add_structure_radio' )
        self.tool_rotate_dihedral_radio = self.builder.get_object ( 'tool_rotate_dihedral_radio' )
        self.element_C_radio   = self.builder.get_object ( 'element_C_radio' )
        self.element_N_radio   = self.builder.get_object ( 'element_N_radio' )
        self.element_O_radio   = self.builder.get_object ( 'element_O_radio' )
        self.element_H_radio   = self.builder.get_object ( 'element_H_radio' )
        # [EN] 2026-09-24, user's own explicit request: "Adicionei os
        # elmentos: B, P, S, F, Cl, Br, I" -- 7 new quick-pick element
        # radios alongside the original C/N/O/H/Other... set (see
        # builder_sidebar.glade's own element_grid, 4 columns x 3 rows).
        self.element_B_radio   = self.builder.get_object ( 'element_B_radio' )
        self.element_P_radio   = self.builder.get_object ( 'element_P_radio' )
        self.element_S_radio   = self.builder.get_object ( 'element_S_radio' )
        self.element_F_radio   = self.builder.get_object ( 'element_F_radio' )
        self.element_Cl_radio  = self.builder.get_object ( 'element_Cl_radio' )
        self.element_Br_radio  = self.builder.get_object ( 'element_Br_radio' )
        self.element_I_radio   = self.builder.get_object ( 'element_I_radio' )
        self.element_other_radio = self.builder.get_object ( 'element_other_radio' )
        self._other_element_symbol = None
        # [EN] 2026-09-24/25: quick-pick fragment row (user's own hand-
        # edited glade -- see _QUICK_FRAGMENT_FILES above). Fetching refs
        # for ALL 12 current fragment_grid buttons, not just the 7 that
        # already have a real .mol2 mapping -- the rest (frag_coo_radio/
        # frag_con_radio1/frag_coo_radio1/frag_coo_radio2) are kept as
        # inert placeholders per the user's own "podemos linkar os
        # fragmentos aos botoes em si mais tarde".
        self.frag_cooh_radio   = self.builder.get_object ( 'frag_cooh_radio' )
        self.frag_coo_radio    = self.builder.get_object ( 'frag_coo_radio' )
        self.frag_con_radio1   = self.builder.get_object ( 'frag_con_radio1' )
        self.frag_nco_radio    = self.builder.get_object ( 'frag_nco_radio' )
        self.frag_ome_radio    = self.builder.get_object ( 'frag_ome_radio' )
        self.frag_hex_radio    = self.builder.get_object ( 'frag_hex_radio' )
        self.frag_pent_radio   = self.builder.get_object ( 'frag_pent_radio' )
        self.frag_benz_radio   = self.builder.get_object ( 'frag_benz_radio' )
        self.frag_furo_radio   = self.builder.get_object ( 'frag_furo_radio' )
        self.frag_coo_radio1   = self.builder.get_object ( 'frag_coo_radio1' )
        self.frag_coo_radio2   = self.builder.get_object ( 'frag_coo_radio2' )
        self.frag_other_radio  = self.builder.get_object ( 'frag_other_radio' )
        self.bond_order_single_radio   = self.builder.get_object ( 'bond_order_single_radio' )
        self.bond_order_double_radio   = self.builder.get_object ( 'bond_order_double_radio' )
        self.bond_order_triple_radio   = self.builder.get_object ( 'bond_order_triple_radio' )
        self.transform_selection_button = self.builder.get_object ( 'transform_selection_button' )
        self.undo_button       = self.builder.get_object ( 'undo_button' )
        self.clean_up_button   = self.builder.get_object ( 'clean_up_button' )
        self.clean_up_adjust_hydrogens_checkbutton = self.builder.get_object ( 'clean_up_adjust_hydrogens_checkbutton' )
        self.optimize_dyff_button = self.builder.get_object ( 'optimize_dyff_button' )

        if getattr ( self.vm_session, "builder_target_object", None ) is None:
            vismol_object = create_empty_vismol_object ( self.vm_session, name = "builder_molecule" )
            self.vm_session.builder_target_object = vismol_object

        # [EN] 2026-09-25, user's own explicit request: "quando um novo
        # vobject e criado para edicao, temos que desativar a
        # visualizacao de todos os outros objetos existentes" -- hides
        # (not discards -- see hide_other_vobjects_for_builder()'s own
        # docstring) every OTHER currently-visible object, so the 3D view
        # stays focused on just what's being edited (this is ALSO what
        # keeps an "Edit in Builder" session's own ORIGINAL vobject from
        # rendering doubled-up with the temp clone now being edited).
        # Covers all 3 real callers of open_window() uniformly (a blank
        # "New" canvas, and both "Edit in Builder" entry points, which
        # set builder_target_object before calling this). Reversed in
        # close_window() below.
        hide_other_vobjects_for_builder ( self.vm_session, self.vm_session.builder_target_object )

        self.vm_session.builder_atom_mode          = True
        self.vm_session.builder_tool               = "add"
        self.vm_session.builder_atom_symbol        = "C"
        self.vm_session.builder_bond_order         = 1
        self.vm_session.builder_bond_aromatic      = False
        self.vm_session.builder_bond_pick_first_atom = None
        self.vm_session.builder_dihedral_pick_atoms  = [ ]
        self.vm_session.builder_dihedral_axis        = None
        self.vm_session.builder_dihedral_subgroup    = None
        self.vm_session.builder_atom_types_labeled_atoms = [ ]
        # [EN] 2026-09-25, user's own explicit request: "quando ele
        # estiver ativo, o novo atomo e o atomo ao qual este novo atomo
        # estiver ligado devem ter seus numeros de hidrogenios
        # atualizados para uma protonacao padrao. Quando desligado, ele
        # deve permitir que, por exemplo, um nitrogenio faca 4 ligacoes
        # com hidrogenios (ions amonio)" -- this checkbox used to be
        # read ONLY at Clean Up click-time (a local widget read, see on_
        # clean_up_button_clicked() below, UNCHANGED). Mirrored here into
        # a live vm_session flag too, so click_mode.py's/atom_ops.py's
        # own AUTOMATIC post-edit hydrogen adjustment (every normal add/
        # bond/replace-element/change-bond-order/attach-fragment
        # operation) can respect the SAME switch -- when OFF, those
        # operations keep whatever bonds the user manually creates
        # exactly as-is, instead of always auto-correcting back to
        # standard valence, letting non-standard states like a 4-bond
        # ammonium nitrogen be built by hand, one bond/click at a time.
        self.vm_session.builder_adjust_hydrogen_count = True
        # [EN] 2026-09-25: elements and fragments now share ONE radio
        # group (element_C_radio's -- see builder_sidebar.glade's
        # element_grid/fragment_grid, user's own explicit request: "o
        # botao 'Add' fica atrelado aos elementos e fragmentos, portanto,
        # apenas um botao dentre os elementos e fragmentos pode estar
        # ativo"), so a fresh Builder session starts on plain Carbon
        # (element_C_radio, the group's own documented default in the
        # glade), NOT a pre-armed fragment -- no fragment is selected
        # until the user actually picks one.
        self.vm_session.builder_selected_fragment = None

        # reflete o estado inicial nos widgets sem disparar os handlers
        # de "toggled" de novo (eles ja rodariam com os valores certos,
        # mas sinalizar aqui deixa explicito o que esta acontecendo)
        self.placemode_toggle.set_active ( True )
        self.placemode_toggle.set_label ( "Editing: ON" )
        self.tool_add_radio.set_active ( True )
        self.element_C_radio.set_active ( True )
        self.bond_order_single_radio.set_active ( True )
        self.clean_up_adjust_hydrogens_checkbutton.set_active ( True )

        self.window.show_all ( )
        self.visible = True

    def close_window ( self, *args ):
        """ Closes the sidebar and turns Builder editing mode back OFF.
        Accepts *args so it can be connected directly to BOTH the
        "Close" button's `clicked` signal (button) -> and the window's
        own `delete-event` (widget, event) -- see builder_sidebar.glade,
        both point at on_close_button_clicked() below, which just calls
        this. Returns True so, when triggered via delete-event (the
        window's own X button), GTK doesn't ALSO try to run its own
        default destroy handling on top of the explicit self.window.
        destroy() already done here.

        [EN] Also discards the Builder's target object if no structural
        edit ever happened during this session (see
        discard_builder_object_if_unused()'s own docstring).
        vm_session.builder_target_object is ALWAYS reset to None below,
        regardless of which case applied -- so the *next* open_window()
        always creates a genuinely fresh blank object, never silently
        "resumes" whatever this session happened to leave behind.

        [EN] 2026-09-25 BUG FIX (user's own report: "a opcao de criar um
        novo sistema pelo builder... so funciona uma vez"): this used to
        leave builder_target_object SET (not None) after a blank-canvas
        session that DID receive a real edit (is_builder_only flips to
        False the moment sync_pdynamo_system() first runs -- see empty_
        object.py) -- neither branch below matched that case, so the
        reset at the bottom never ran. The NEXT time the user asked for
        a brand-new blank system (builder_entry_dialog.choose_and_open_
        builder()'s "New" option), open_window()'s own "only create a
        blank object when builder_target_object is still None" guard
        would find it very much NOT None -- and silently reopen the
        PREVIOUS (already-committed, already-promoted) system instead of
        starting fresh, matching the user's exact "works once" report.
        Checked every caller of open_window() before making this
        unconditional: the ONLY place that legitimately wants a "resume"
        (empty_object.begin_editing_existing_system()) always explicitly
        sets builder_target_object itself, immediately before calling
        open_window() -- completely independent of whatever this method
        leaves behind -- so there is no real caller left to preserve the
        old "leave it set" behaviour for. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is not None and getattr ( target_object, "is_builder_only", False ):
            discard_builder_object_if_unused ( target_object )
        elif target_object is not None and getattr ( target_object, "builder_edit_source_e_id", None ) is not None:
            # [EN] "Edit in Builder" session on an already-loaded system
            # (see empty_object.begin_editing_existing_system()) -- unlike
            # the blank-canvas case above, this object already has a real,
            # committed e_id (from clone_system()), so there's no "was it
            # ever promoted" question; finish_editing_existing_system()
            # decides whether to fold it back onto the original (atom
            # count unchanged) or keep it as its own new system (atom
            # count changed).
            finish_editing_existing_system ( target_object )
        self.vm_session.builder_target_object = None

        # [EN] 2026-09-25: counterpart to hide_other_vobjects_for_
        # builder() in open_window() above -- makes every object this
        # session hid visible again (whatever happened to target_object
        # itself above is unrelated to this; some OTHER vobject the user
        # was simply not editing may have been hidden purely so the 3D
        # view stayed focused during this session).
        restore_hidden_vobjects_for_builder ( self.vm_session )

        self.vm_session.builder_atom_mode = False
        self.vm_session.builder_bond_pick_first_atom = None
        self.vm_session.builder_dihedral_pick_atoms  = [ ]
        self.vm_session.builder_dihedral_axis        = None
        self.vm_session.builder_dihedral_subgroup    = None
        # [EN] 2026-09-24: the Dihedral Angle window is now its own
        # top-level window (see dihedral_angle_window.py), not a widget
        # nested inside this one -- unlike the old inline slider (which
        # got destroyed for free along with this window), it needs an
        # explicit close here, or it would be left dangling open after
        # the Builder itself closes.
        self.disarm_dihedral_slider ( )
        # [EN] The Atom Types window's own 3D-view index-label overlay
        # (see atom_types_window.py's _refresh_rows()) points at atoms
        # belonging to the object THIS Builder session was editing --
        # stale/dangling once that session ends, regardless of whether
        # the Atom Types window itself is still open or already closed.
        self.vm_session.builder_atom_types_labeled_atoms = [ ]
        if getattr ( self, "window", None ) is not None:
            self.window.destroy ( )
        self.visible = False
        return True

    # ------------------------------------------------------------------
    #  Signal handlers (referenced by name in builder_sidebar.glade)
    # ------------------------------------------------------------------

    def on_placemode_toggled ( self, button ):
        """ The lighter-weight PAUSE toggle -- see this module's design
        note at the top for how this differs from actually closing the
        window. """
        is_active = button.get_active ( )
        self.vm_session.builder_atom_mode = is_active
        button.set_label ( "Editing: ON" if is_active else "Editing: OFF" )

    def on_tool_changed ( self, button ):
        """ [EN] GtkRadioButton fires "toggled" for BOTH the button that
        just became active AND the one that just became inactive -- only
        act on the one reporting active=True, otherwise this runs TWICE
        per click, the second time with the (now wrong) previous tool.

        [EN] 2026-09-25: there is no longer a dedicated "Fragment" TOOL
        radio (user's own explicit request -- "o botao 'Add' fica
        atrelado aos elementos e fragmentos"). vm_session.builder_tool
        can still become "attach_fragment" internally (that's the exact
        string value vismol_glcore.py's render() hook dispatches on --
        see handle_click_to_attach_fragment()), it's just set by on_
        element_changed()/on_fragment_quick_changed() below (whichever
        of the merged element/fragment radio group was picked last)
        instead of by a click on THIS tool_add_radio branch. Clicking
        "Add" directly here always resets to plain "add" (single-atom)
        mode, even if a fragment still looks selected in the row below --
        a deliberate, simple manual override, not yet unified further. """
        if not button.get_active ( ):
            return
        if button is self.tool_add_radio:
            self.vm_session.builder_tool = "add"
        elif button is self.tool_delete_radio:
            self.vm_session.builder_tool = "delete"
        elif button is self.tool_move_radio:
            self.vm_session.builder_tool = "move"
        elif button is self.tool_bond_order_radio:
            self.vm_session.builder_tool = "bond_order"
        elif button is self.tool_add_structure_radio:
            self.vm_session.builder_tool = "add_structure"
            self.main.structure_library_window.open_window ( )
        elif button is self.tool_rotate_dihedral_radio:
            self.vm_session.builder_tool = "rotate_dihedral"

        # [EN] Leaving the "Rotate Dihedral" tool (or never having been on
        # it) drops any in-progress 4-atom pick and disarms the slider --
        # same stale-state reasoning as the Bond Order block below, just
        # for this tool's own armed state.
        if self.vm_session.builder_tool != "rotate_dihedral":
            self.vm_session.builder_dihedral_pick_atoms = [ ]
            self.vm_session.builder_dihedral_axis       = None
            self.vm_session.builder_dihedral_subgroup   = None
            self.disarm_dihedral_slider ( )

        # [EN] Leaving the "Bond Order" tool (or never having been on it)
        # drops any pending "first atom armed, waiting for its partner"
        # click_mode.handle_click_to_set_bond_order() may have left behind
        # -- otherwise switching to Add/Delete/Move mid-pick and back
        # could silently resume bonding to an atom the user has long since
        # forgotten they clicked (and, worse, one that clean-up/undo might
        # have since removed -- handle_click_to_set_bond_order() already
        # guards against a STALE reference, but not against a stale
        # SURVIVING one the user simply doesn't remember picking).
        if self.vm_session.builder_tool != "bond_order":
            self.vm_session.builder_bond_pick_first_atom = None

    def on_element_changed ( self, button ):
        """ Same "only act on the newly-active one" reasoning as
        on_tool_changed() above.

        [EN] "Other..." (user's own request -- a button to pick ANY
        element, not just the 4 quick-pick ones here) opens the full
        periodic table picker (periodic_table_dialog.choose_element())
        the moment it becomes the active radio, rather than sitting
        there as its own separate, silent "mode" -- picking a symbol
        sets builder_atom_symbol AND relabels this radio to show what
        was picked (so it stays visible which element is actually
        armed); cancelling with nothing EVER picked yet on this radio
        falls back to Carbon (the same default the sidebar opens with)
        instead of leaving an ambiguous "Other..." selected with no
        real symbol behind it -- cancelling a RE-pick (something was
        already chosen here before) just keeps that previous choice
        active, unchanged.

        [EN] 2026-09-25, user's own explicit request: "o botao 'Add'
        fica atrelado aos elementos e fragmentos" -- element_grid now
        shares ONE radio group with fragment_grid (see builder_sidebar.
        glade), so picking an element here also means "Add" is the
        governing tool again (overriding whatever "attach_fragment" tool
        state a previously-picked fragment may have left armed). Synced
        FIRST, before the final builder_tool assignment below, since
        tool_add_radio.set_active() can itself fire on_tool_changed()
        (which sets builder_tool = "add") -- setting builder_tool
        explicitly afterwards keeps this method's own value as the one
        that actually sticks, regardless of signal-firing order. """
        if not button.get_active ( ):
            return

        if button is self.element_other_radio:
            from gui.windows.builder.periodic_table_dialog import choose_element
            symbol = choose_element ( self.main )
            if symbol is not None:
                self._other_element_symbol = symbol
                self.element_other_radio.set_label ( symbol )
                self.tool_add_radio.set_active ( True )
                self.vm_session.builder_atom_symbol = symbol
                self.vm_session.builder_tool = "add"
            elif self._other_element_symbol is None:
                self.element_C_radio.set_active ( True )   # nothing chosen yet -- fall back to the default
            return

        symbol_by_button = {
            self.element_C_radio:  "C",
            self.element_N_radio:  "N",
            self.element_O_radio:  "O",
            self.element_H_radio:  "H",
            self.element_B_radio:  "B",
            self.element_P_radio:  "P",
            self.element_S_radio:  "S",
            self.element_F_radio:  "F",
            self.element_Cl_radio: "Cl",
            self.element_Br_radio: "Br",
            self.element_I_radio:  "I",
        }
        symbol = symbol_by_button.get ( button )
        if symbol is not None:
            self.tool_add_radio.set_active ( True )
            self.vm_session.builder_atom_symbol = symbol
            self.vm_session.builder_tool = "add"

    def on_fragment_quick_changed ( self, button ):
        """ [EN] 2026-09-24/25, user's own explicit request: "vamos
        colocar abaixo do elementos os fragmentos mais importantes na
        forma de botoes (com isso vamos deixar janela de fragmentos
        apenas para quando o usuario chamar o botao 'other')" -- and,
        2026-09-25: "o botao 'Add' fica atrelado aos elementos e
        fragmentos, portanto, apenas um botao dentre os elementos e
        fragmentos pode estar ativo" -- fragment_grid now shares ONE
        radio group with element_grid (see builder_sidebar.glade), so
        picking a fragment here automatically deactivates whichever
        element radio was active, and vice versa -- pure native GTK
        group exclusivity, no extra bookkeeping needed for that part.
        Same "only act on the newly-active one" reasoning as on_tool_
        changed()/on_element_changed() above.

        Picking one of the mapped quick fragments (see this module's own
        _QUICK_FRAGMENT_FILES) mirrors exactly what fragment_library_
        window.py's own _load_and_select() already does when a fragment
        is chosen from the full Library (set builder_selected_fragment,
        force builder_tool to "attach_fragment", and sync the Tool
        grid's own "Add" radio -- there is no separate "Fragment" tool
        radio anymore) -- so a quick-pick click is immediately ready to
        attach, no extra step. Buttons with no entry in _QUICK_FRAGMENT_
        FILES (the user's own newer, not-yet-linked ones -- see that
        dict's own comment) simply no-op here, same as clicking nothing.
        "Other..." does none of the fragment-loading itself -- it only
        opens the full Fragment Library window, which sets all of the
        above the moment something is actually chosen inside it. """
        if not button.get_active ( ):
            return

        if button is self.frag_other_radio:
            self.main.fragment_library_window.open_window ( )
            return

        filename = None
        for attr_name, mol2_name in _QUICK_FRAGMENT_FILES.items ( ):
            if button is getattr ( self, attr_name, None ):
                filename = mol2_name
                break
        if filename is None:
            return

        file_path = os.path.join ( self.main.home, "src/gui/windows/builder/fragments/functional_groups", filename )
        try:
            fragment = load_fragment ( file_path )
        except FragmentError as error:
            self.main.simple_dialog.error ( msg = str ( error ) )
            return

        # [EN] Sync "Add" active FIRST -- it may itself fire on_tool_
        # changed() (which sets builder_tool = "add") -- so the more
        # specific "attach_fragment" value set below always wins as the
        # final state, regardless of signal-firing order.
        self.tool_add_radio.set_active ( True )
        self.vm_session.builder_selected_fragment = fragment
        self.vm_session.builder_tool = "attach_fragment"
        if getattr ( self.main, "statusbar_main", None ) is not None:
            self.main.statusbar_main.push ( 1, "Selected: {} -- click a hydrogen atom in the Builder to attach it.".format ( fragment["name"] ) )

    def on_bond_order_changed ( self, button ):
        """ Same "only act on the newly-active one" reasoning as
        on_tool_changed()/on_element_changed() above. This selection
        applies to bonds CREATED from now on (drag-to-bond, 'b' key), to
        an EXISTING bond via Ctrl+click on it (click_mode.apply_
        selected_bond_order()), and to an existing bond via the "Bond
        Order" click-tool above (click_mode.handle_click_to_set_bond_
        order()).

        [EN] 2026-09-25, user's own explicit request: "vamos descartar
        as ligacoes tipo aromaticas" -- the "Aromatic" option is gone
        from this row (bond_order_grid now only has Single/Double/
        Triple). builder_bond_aromatic stays False from here on (its
        init in open_window() is unchanged) -- other code that still
        reads it (click_mode.py/vismol_glcore.py, e.g. drag-distance
        preview colouring, aromatic PERCEPTION on import) all does so
        via getattr(..., False), so leaving the variable itself in place
        but permanently False from the UI's side is safe, not a partial
        removal. """
        if not button.get_active ( ):
            return
        if button is self.bond_order_single_radio:
            self.vm_session.builder_bond_order    = 1
            self.vm_session.builder_bond_aromatic = False
        elif button is self.bond_order_double_radio:
            self.vm_session.builder_bond_order    = 2
            self.vm_session.builder_bond_aromatic = False
        elif button is self.bond_order_triple_radio:
            self.vm_session.builder_bond_order    = 3
            self.vm_session.builder_bond_aromatic = False

    def on_adjust_hydrogen_count_toggled ( self, button ):
        """ [EN] 2026-09-25, user's own explicit request -- see open_
        window()'s own comment on vm_session.builder_adjust_hydrogen_
        count for the full reasoning. Keeps that live flag in sync with
        the checkbox's own state, read directly by click_mode.py's/
        atom_ops.py's automatic post-edit hydrogen adjustment (on_clean_
        up_button_clicked() below still reads the widget directly at
        click-time, unchanged -- both always agree since it's the SAME
        checkbox). """
        self.vm_session.builder_adjust_hydrogen_count = button.get_active ( )

    def on_atom_types_button_clicked ( self, button ):
        """ Opens the Atom Types window (atom_types_window.
        AtomTypesWindow, instantiated once on main_window -- same
        pattern as fragment_library_window above) on whichever atoms are
        currently in the app's own VIEWING selection (shift+click/
        shift+drag -- a separate, pre-existing mechanism from any
        Builder click tool, confirmed to already work alongside them
        with no conflict) -- filtered to atoms belonging to THIS
        Builder's own target object, since a viewing selection isn't
        scoped to any one object and could include atoms from something
        else entirely on screen. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return

        selection = self.vm_session.selections.get ( self.vm_session.current_selection )
        selected_atoms = getattr ( selection, "selected_atoms", None ) or set ( )
        atom_ids = { atom.atom_id for atom in selected_atoms if atom.vm_object is target_object }

        if not atom_ids:
            if getattr ( self.main, "statusbar_main", None ) is not None:
                self.main.statusbar_main.push ( 1, "Atom Types: select some atoms first (shift+click or shift+drag)." )
            return

        self.main.atom_types_window.open_window ( target_object, atom_ids )

    def on_transform_selection_button_clicked ( self, button ):
        """ Opens the Transform Selection window (transform_selection_
        window.TransformSelectionWindow, instantiated once on main_window)
        on whichever atoms are currently in the app's own VIEWING
        selection -- same filtering-to-the-Builder's-own-target-object
        reasoning as on_atom_types_button_clicked() just above. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return

        selection = self.vm_session.selections.get ( self.vm_session.current_selection )
        selected_atoms = getattr ( selection, "selected_atoms", None ) or set ( )
        atom_ids = { atom.atom_id for atom in selected_atoms if atom.vm_object is target_object }

        if not atom_ids:
            if getattr ( self.main, "statusbar_main", None ) is not None:
                self.main.statusbar_main.push ( 1, "Transform Selection: select some atoms first (shift+click or shift+drag)." )
            return

        self.main.transform_selection_window.open_window ( target_object, atom_ids )

    def arm_dihedral_slider ( self, angle_deg ):
        """ [EN] 2026-09-24: the inline slider this method used to drive
        directly was extracted into its own window (dihedral_angle_
        window.DihedralAngleWindow -- user's own explicit request, "a
        edicao do diedro deve ser feita em uma janela a parte, como
        acontece em 'Transform Selections'") -- this now just opens/
        updates that window instead. Still called the exact same way by
        click_mode.handle_click_to_pick_dihedral_atom() (via vm_session.
        main.builder_sidebar_window, the same cross-window-sync pattern
        fragment_library_window.py already uses for THIS sidebar's Tool
        radio), so no change needed at that call site. """
        self.main.dihedral_angle_window.open_window ( angle_deg )

    def disarm_dihedral_slider ( self ):
        """ Counterpart to arm_dihedral_slider() above -- closes the
        Dihedral Angle window (tool switched away, pick cancelled, or
        this sidebar itself closing). A no-op if it's already closed
        (e.g. the user closed it by hand mid-edit). """
        dihedral_window = getattr ( self.main, "dihedral_angle_window", None )
        if dihedral_window is not None and dihedral_window.visible:
            dihedral_window.close_window ( )

    def on_undo_button_clicked ( self, button ):
        """ Calls atom_ops.undo() on the CURRENT target object -- see
        that function's own docstring (and push_undo_snapshot()'s) for
        exactly which actions this can step back through (place/replace
        an atom, a whole drag-to-create-a-bonded-atom gesture, deleting
        an atom or bond, changing a bond's order). """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return
        atom_ops_undo ( target_object )

    def on_clean_up_button_clicked ( self, button ):
        """ Runs atom_ops.clean_up_structure() on the WHOLE current
        target object (atom_ids=None -- see that function's own
        docstring for the difference between whole-molecule and
        localised modes; this button always uses whole-molecule, per
        explicit request). Pushes an undo snapshot first, same as every
        other Builder mutation, since this can move every atom in the
        object at once.

        adjust_hydrogen_count comes straight from the sidebar's own
        checkbox (user's own request) -- when unticked, Clean Up only
        relaxes heavy-atom positions/angles/planarity and leaves
        whichever hydrogens already exist exactly where they are,
        neither repositioned nor added/removed (see clean_up_structure()'s
        own docstring for the full reasoning). """
        from gui.windows.builder.atom_ops import push_undo_snapshot
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return
        adjust_hydrogen_count = ( self.clean_up_adjust_hydrogens_checkbutton is None
                                   or self.clean_up_adjust_hydrogens_checkbutton.get_active ( ) )
        push_undo_snapshot ( target_object )
        clean_up_structure ( target_object, atom_ids = None, adjust_hydrogen_count = adjust_hydrogen_count )

        from gui.windows.builder.empty_object import sync_pdynamo_system
        sync_pdynamo_system ( target_object )

    def on_optimize_dyff_button_clicked ( self, button ):
        """ User's own request: "otimizar a geometria com o DYFF, vamos
        criar um sistema provisorio no background, otimizar por 100 passos
        e atualizar as coordenadas, depois podemos descartar o sistema
        provisorio" -- runs atom_ops.optimize_geometry_dyff() (a REAL
        force-field minimisation on a disposable scratch pDynamo system,
        see that function's own docstring) on the WHOLE current target
        object. Pushes an undo snapshot first, same as Clean Up, since
        this can move every atom in the object at once. """
        from gui.windows.builder.atom_ops import push_undo_snapshot
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return
        push_undo_snapshot ( target_object )
        ok, message = optimize_geometry_dyff ( target_object )
        if ok:
            self.main.statusbar_main.push ( 1, message )
            from gui.windows.builder.empty_object import sync_pdynamo_system
            sync_pdynamo_system ( target_object )
        else:
            self.main.simple_dialog.error ( msg = message )

    def on_close_button_clicked ( self, *args ):
        return self.close_window ( )
