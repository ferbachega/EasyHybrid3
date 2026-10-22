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
#
#  [EN] 2026-10-17, user's request (LibreOffice-like layout): there is no
#  separate window any more. The tools are in the main window's
#  toolbar_edit (MainWindow.glade), shown only while editing, and the
#  options of the active tool (elements + fragments for "Add", bond order
#  for "Bond Order") in the right side panel (MainWindow.show_sidebar_page()).
#  builder_sidebar.glade is no longer loaded; the class keeps its name and
#  its open_window()/close_window()/visible interface for its callers.
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
from gui.windows.builder.atom_ops     import redo as atom_ops_redo
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
# [EN] 2026-10-02, user's own glade edit: every fragment button now has its
# meaning in its TOOLTIP (the widget ids are older names and were kept as
# the user left them -- e.g. frag_pent_radio is now Nitro, frag_hex_radio
# is Cyclopentadiene; they used to load pentyl/hexyl). Paths are relative
# to fragments/ (category/file). Formal charges (nitro, carboxylate) and the
# hypervalent S/P of sulfone/phosphonic acid are supported by
# fragment_library.load_fragment() since this same change.
_QUICK_FRAGMENT_FILES = {
    "frag_cooh_radio":         "functional_groups/carboxyl.mol2",        # Carboxylic acid, -COOH
    "frag_coo_radio":          "functional_groups/carboxylate.mol2",     # Carboxylate, -COO(-)
    "frag_nco_radio":          "functional_groups/formamido.mol2",       # Amide, -N-C=O (N-linked)
    "frag_con_radio1":         "functional_groups/amide.mol2",           # Primary amide, -C(=O)NH2
    "frag_coo_radio1":         "functional_groups/methylsulfonyl.mol2",  # Sulfone, -SO2(-CH3)
    "frag_coo_radio2":         "functional_groups/phosphono.mol2",       # -PO3 (phosphonic acid, -PO3H2)
    "frag_ome_radio":          "functional_groups/methoxy.mol2",         # Methoxy, -O-CH3
    "frag_pent_radio":         "functional_groups/nitro.mol2",           # Nitro, -NO2
    "frag_benz_radio":         "functional_groups/phenyl.mol2",          # Phenyl ring (Kekulized)
    "frag_chexane_radio":      "rings/cyclohexyl.mol2",                  # Cyclohexane
    "frag_hex_radio":          "rings/cyclopentadienyl.mol2",            # Cyclopentadiene
    # tooltip says "Cyclopropane" (same as frag_cpropane_radio3), but its
    # icon is frag_cpent.png -- mapped to cyclopentyl, see the note to the user.
    "frag_ciclopropane_radio": "rings/cyclopentyl.mol2",
    "frag_cbutane_radio2":     "rings/cyclobutyl.mol2",                  # Cyclobutane
    "frag_cpropane_radio3":    "rings/cyclopropyl.mol2",                 # Cyclopropane
}


# [EN] 2026-10-17: Force Field Settings (toolbar_edit "FF" toggle, side
# panel page builder_ff_settings), kept in .config.json as
# gl_parameters['builder_ff_settings']. Defaults = the behaviour before
# the settings existed (ask on close, plain DYFF).
FF_SETTINGS_DEFAULTS = {
    "force_field": "dyff",        # none | dyff | dyff_xtb | opls
    "method":      "gfn2",        # gfn2 | gxtb (dyff_xtb only)
    "water":       True,          # implicit water in the xTB calculations
    "bonded":      True,          # Hessian bonds/angles + torsions
    "cm5_scale":   1.0,
    "nb_model":    "cutoff",      # cutoff | full
    "on_close":    "ask",         # ask | apply
    "opls_cm5_scale": 1.0,        # GFN1-xTB CM5 scale of the OPLS ligand charges (opls only)
}

FF_NAMES = { "none": "no force field", "dyff": "DYFF", "dyff_xtb": "DYFF + xTB refinement", "opls": "OPLS-AA" }

FF_DESCRIPTIONS = {
    "none":     "The system is kept without a force field.",
    "dyff":     "Atom types from element + bond orders. No partial charges, except those set in DYFF Parameters.",
    "dyff_xtb": "DYFF with CM5 charges, Hessian bonds/angles and torsions from xTB (standard residues from the "
                "residue library). Runs in the background after the Builder closes.",
    "opls":     "OPLS-AA parameter set. If parameters are missing, a dialog offers the Prepare OPLS System procedure: "
                "OPLS-AA classes by rules (Lennard-Jones, bonded terms), CM5 charges from one GFN1-xTB calculation, "
                "missing terms estimated (Dodda et al., JCTC 2015). No solvation: that is a separate step.",
}


class BuilderSidebarWindow ( ):
    """ Builder editing sidebar. """

    def __init__ ( self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.visible    = False
        self._widgets_bound = False

    def open_window ( self ):
        """ Opens the Builder -- and, per this module's own design note
        above, turns Builder editing mode ON, creating a fresh empty
        target object first if none exists yet. [EN] 2026-10-17: no
        separate window any more -- shows the main window's toolbar_edit
        and the right side panel (see _bind_widgets()). """
        if self.visible:
            return

        self._bind_widgets ( )

        if getattr ( self.vm_session, "builder_target_object", None ) is None:
            vismol_object = create_empty_vismol_object ( self.vm_session, name = "builder_molecule" )
            self.vm_session.builder_target_object = vismol_object
        self._start_session ( )

    def _bind_widgets ( self ):
        """ [EN] 2026-10-17, user's request (LibreOffice-like layout): the
        Builder's tools live in the main window -- MainWindow.glade's
        toolbar_edit (visible only while editing) and, for the tool that
        needs it, a page of the right side panel (sidebar_root_box /
        sidebar_stack: elements + fragments for "Add", Single/Double/Triple
        for "Bond Order"). Same widget ids as the old builder_sidebar.glade
        window (no longer loaded); the main window's glade has no signals for
        them, so they are connected here, once. """
        if self._widgets_bound:
            return
        self.builder = self.main.builder
        self.toolbar_edit = self.builder.get_object ( 'toolbar_edit' )

        self.placemode_toggle  = self.builder.get_object ( 'placemode_toggle' )
        self.tool_add_radio    = self.builder.get_object ( 'tool_add_radio' )
        self.tool_delete_radio = self.builder.get_object ( 'tool_delete_radio' )
        self.tool_move_radio   = self.builder.get_object ( 'tool_move_radio' )
        self.tool_bond_order_radio = self.builder.get_object ( 'tool_bond_order_radio' )
        # [EN] 2026-09-25: no more dedicated "Fragment" TOOL radio -- see
        # this module's own _QUICK_FRAGMENT_FILES comment above.
        # [EN] 2026-10-17: no "Structure" tool radio either -- structures are
        # the "Structures" tab of the Add panel (see _bind_add_panel_tabs()).
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
        # [EN] 2026-10-02: every mapped quick-fragment button gets its
        # widget attribute and its "toggled" handler from here (the Python
        # connections below are the only ones -- 2026-10-17).
        for attr_name in _QUICK_FRAGMENT_FILES:
            widget = self.builder.get_object ( attr_name )
            setattr ( self, attr_name, widget )
            if widget is not None:
                widget.connect ( "toggled", self.on_fragment_quick_changed )
        self._last_quick_fragment = ( None, None )
        self.bond_order_single_radio   = self.builder.get_object ( 'bond_order_single_radio' )
        self.bond_order_double_radio   = self.builder.get_object ( 'bond_order_double_radio' )
        self.bond_order_triple_radio   = self.builder.get_object ( 'bond_order_triple_radio' )
        self.transform_selection_button = self.builder.get_object ( 'transform_selection_button' )
        self.undo_button       = self.builder.get_object ( 'undo_button' )
        self.redo_button       = self.builder.get_object ( 'redo_button' )
        # [EN] 2026-10-03 user request: the Bond Order radios only serve the
        # "Bond Order" tool (re-defining an EXISTING bond), so their row is
        # shown only while that tool is active (see _sync_side_panel()).
        # "Add" no longer reads them: a dragged bond's order always comes
        # from the drag distance (click_mode.update/finish_bond_drag()).
        # [EN] 2026-10-03: the Clean Up button was removed (user request);
        # atom_ops.clean_up_structure() itself is kept.
        # [EN] 2026-10-02 user request ("ferramenta para adicionar os
        # hidrogenios corretamente... incorporados aos residuos"): an "Add
        # Hydrogens" button right below Optimize. [EN] 2026-10-03: now
        # defined in builder_sidebar.glade (id add_hydrogens_button, handler
        # on_add_hydrogens_button_clicked); _create_add_hydrogens_button()
        # just returns it, creating one only if an older glade lacks it.
        self.optimize_dyff_button = self.builder.get_object ( 'optimize_dyff_button' )
        self.add_hydrogens_button = self._create_add_hydrogens_button ( )
        # dark interface theme (src/gui/theme.py): lighten the black line-
        # drawing icons (ring fragments...) so they stay visible; no-op with
        # the default theme
        try:
            from gui.theme import adapt_icons_for_theme
            adapt_icons_for_theme ( self.builder.get_object ( 'sidebar_root_box' ) )
        except Exception:
            pass
        self.clean_up_adjust_hydrogens_checkbutton = self.builder.get_object ( 'clean_up_adjust_hydrogens_checkbutton' )

        handlers = {
            'placemode_toggle':                      ( 'toggled', self.on_placemode_toggled ),
            'clean_up_adjust_hydrogens_checkbutton': ( 'toggled', self.on_adjust_hydrogen_count_toggled ),
            'undo_button':                           ( 'clicked', self.on_undo_button_clicked ),
            'redo_button':                           ( 'clicked', self.on_redo_button_clicked ),
            'optimize_dyff_button':                  ( 'clicked', self.on_optimize_dyff_button_clicked ),
            'add_hydrogens_button':                  ( 'clicked', self.on_add_hydrogens_button_clicked ),
            'atom_types_button':                     ( 'clicked', self.on_atom_types_button_clicked ),
            'solvate_button':                        ( 'clicked', self.on_solvate_button_clicked ),
            'close_button':                          ( 'clicked', self.on_close_button_clicked ),
            'frag_other_radio':                      ( 'toggled', self.on_fragment_quick_changed ),
        }
        for name in ( 'tool_add_radio', 'tool_delete_radio', 'tool_bond_order_radio', 'tool_move_radio',
                      'tool_rotate_dihedral_radio', 'transform_selection_button' ):
            handlers[name] = ( 'toggled', self.on_tool_changed )
        for name in ( 'bond_order_single_radio', 'bond_order_double_radio', 'bond_order_triple_radio' ):
            handlers[name] = ( 'toggled', self.on_bond_order_changed )
        for symbol in ( 'C', 'N', 'O', 'H', 'B', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'other' ):
            handlers['element_{}_radio'.format ( symbol )] = ( 'toggled', self.on_element_changed )
        for name, ( signal, handler ) in handlers.items ( ):
            widget = self.builder.get_object ( name )
            if widget is not None:
                widget.connect ( signal, handler )
        self._bind_add_panel_tabs ( )
        self._bind_delete_panel ( )
        self._bind_dihedral_panel ( )
        self._bind_transform_panel ( )
        self._bind_ff_panel ( )
        self._widgets_bound = True

    # ------------------------------------------------------------------
    #  Side panel pages (MainWindow.glade's sidebar_stack) -- 2026-10-17
    # ------------------------------------------------------------------

    def _bind_add_panel_tabs ( self ):
        """ [EN] 2026-10-17, user's request: the Structure Library is the
        "Structures" tab of the Add panel, next to elements and fragments
        ("Chemical"), instead of a tool of its own with a separate window.
        Picking a row previews it and arms the "add_structure" tool (the
        next click in the view places it -- click_mode.handle_click_to_add_
        structure()); going back to "Chemical" restores the element or
        fragment that was armed there. """
        self.add_notebook            = self.builder.get_object ( 'builder_add_notebook' )
        self.structures_page         = self.builder.get_object ( 'builder_structures_page' )
        self.structure_treestore     = self.builder.get_object ( 'builder_structure_treestore' )
        self.structure_treeview      = self.builder.get_object ( 'builder_structure_treeview' )
        self.structure_preview_box   = self.builder.get_object ( 'builder_structure_preview_box' )
        self.structure_status_label  = self.builder.get_object ( 'builder_structure_status_label' )
        self.structure_preview       = None   # created the first time the tab is shown (GL widget)
        self._selected_structure     = None
        self._chemical_state         = None
        self.add_notebook.connect ( 'switch-page', self.on_add_notebook_switch_page )
        self.structure_treeview.get_selection ( ).connect ( 'changed', self.on_structure_selection_changed )
        self.structure_treeview.connect ( 'row-activated', self.on_structure_row_activated )

    def _populate_structure_tree ( self ):
        """ Re-scanned at every Builder session (a structure added to the
        folder while EasyHybrid is running shows up next time), every
        category collapsed -- same as the old Structure Library window. """
        from gui.windows.builder.structure_library import discover_structure_files
        self.structure_treestore.clear ( )
        for category_name, entries in discover_structure_files ( ):
            category_iter = self.structure_treestore.append ( None, [ category_name, "" ] )
            for display_name, file_path in entries:
                self.structure_treestore.append ( category_iter, [ display_name, file_path ] )
        self.structure_treeview.collapse_all ( )

    def _structures_tab_active ( self ):
        return self.add_notebook.get_current_page ( ) == self.add_notebook.page_num ( self.structures_page )

    def on_add_notebook_switch_page ( self, notebook, page, page_num ):
        """ "switch-page" runs before the current page changes: `page` is the
        new one. """
        if page is self.structures_page:
            if self.structure_preview is None:
                from gui.windows.builder.molecule_preview import MoleculePreview
                self.structure_preview = MoleculePreview ( )
                self.structure_preview_box.pack_start ( self.structure_preview.widget, True, True, 0 )
                self.structure_preview.widget.show ( )
            tool = getattr ( self.vm_session, "builder_tool", None )
            if tool in ( "add", "attach_fragment" ):
                self._chemical_state = ( tool, getattr ( self.vm_session, "builder_atom_symbol", "C" ),
                                         getattr ( self.vm_session, "builder_selected_fragment", None ) )
            if self.tool_add_radio.get_active ( ):
                self.vm_session.builder_tool = "add_structure"
            if self._selected_structure is None:
                self.structure_status_label.set_text ( "No structure selected." )
        elif self._chemical_state is not None and self.tool_add_radio.get_active ( ):
            tool, symbol, fragment = self._chemical_state
            self.vm_session.builder_tool = tool
            self.vm_session.builder_atom_symbol = symbol
            self.vm_session.builder_selected_fragment = fragment

    def on_structure_selection_changed ( self, selection ):
        """ Single click: preview + arm (the tab itself means "add a
        structure", no separate Use step). Category rows and malformed files
        only clear the preview. """
        from gui.windows.builder.structure_library import load_structure, StructureError
        model, tree_iter = selection.get_selected ( )
        file_path = model[tree_iter][1] if tree_iter is not None else ""
        if not file_path:
            if self.structure_preview is not None:
                self.structure_preview.clear ( )
            return
        try:
            structure = load_structure ( file_path )
        except StructureError as error:
            if self.structure_preview is not None:
                self.structure_preview.clear ( )
            self.structure_status_label.set_text ( "Error: {}".format ( error ) )
            return
        if self.structure_preview is not None:
            self.structure_preview.show ( structure["name"], structure["atoms"], structure["bonds"] )
        self._selected_structure = structure
        self.vm_session.builder_selected_structure = structure
        if self.tool_add_radio.get_active ( ):
            self.vm_session.builder_tool = "add_structure"
        # . 2026-10-21: total charge and multiplicity (metal, oxidation state, spin note)
        from gui.windows.builder.structure_library import describe_structure, describe_structure_text
        try:
            info = describe_structure ( structure, file_path )
            structure["total_charge"], structure["multiplicity"] = info["charge"], info["multiplicity"]
            text = describe_structure_text ( structure, file_path )
        except Exception as error:
            text = "{} ({} atoms)".format ( structure["name"], len ( structure["atoms"] ) )
        self.structure_status_label.set_text ( text + "\nClick in the view to place it." )

    def on_structure_row_activated ( self, treeview, path, column ):
        """ Double-click on a category row expands/collapses it. """
        if not self.structure_treestore[path][1]:
            if treeview.row_expanded ( path ):
                treeview.collapse_row ( path )
            else:
                treeview.expand_row ( path, False )

    def _bind_transform_panel ( self ):
        """ [EN] 2026-10-17, user's request: Transform Selection is a tool
        (toolbar radio) with its sliders in the side panel instead of a
        window, and it sets itself up -- entering it pauses Editing (the
        viewing selection only works with Editing OFF, see on_placemode_
        toggled()) and switches the selection to viewing mode; leaving it
        for another tool turns Editing back on. The sliders act on the
        CURRENT viewing selection (atoms of the Builder object), re-read at
        the start of every drag; the header follows it live (see
        _update_undo_redo_sensitivity()). Same rules as the old
        TransformSelectionWindow: live sliders, delta from the last value,
        rotation around the selection centroid, one undo snapshot per drag,
        sliders back to 0 on release. """
        self.transform_header_label = self.builder.get_object ( 'builder_transform_header_label' )
        self.transform_scales = { }
        for key in ( "tx", "ty", "tz", "rx", "ry", "rz" ):
            scale = self.builder.get_object ( 'builder_transform_{}_scale'.format ( key ) )
            self.transform_scales[key] = scale
            scale.connect ( 'value-changed', self.on_transform_slider_changed, key )
            scale.connect ( 'button-release-event', self.on_transform_slider_released )
        self._transform_last_values = { key: 0.0 for key in self.transform_scales }
        self._transform_suppress    = False
        self._transform_gesture_atoms = None
        self._transform_paused_editing = False
        self._transform_sync_timer  = None

    def _bind_ff_panel ( self ):
        """ [EN] 2026-10-17, user's request: "Force Field Settings" button --
        the force field (and its options) for the systems drawn in the
        Builder. A toggle: while on, the side panel shows these settings
        instead of the tool's page ("Done" or picking a tool turns it off).
        Used by close_window() -> _ask_incorporate_with_dyff() /
        _incorporate_with_dyff(). """
        b = self.builder.get_object
        self.ff_settings_toggle = b ( 'ff_settings_toggle' )
        self.ff_radios = { "none": b ( 'ff_none_radio' ), "dyff": b ( 'ff_dyff_radio' ),
                           "dyff_xtb": b ( 'ff_dyff_xtb_radio' ), "opls": b ( 'ff_opls_radio' ) }
        self.ff_widgets = { "method": b ( 'ff_method_combo' ), "water": b ( 'ff_water_check' ),
                            "bonded": b ( 'ff_bonded_check' ), "cm5_scale": b ( 'ff_cm5_spin' ),
                            "nb_model": b ( 'ff_nb_combo' ), "on_close": b ( 'ff_on_close_combo' ),
                            "opls_cm5_scale": b ( 'ff_opls_cm5_spin' ) }
        self._ff_loading = True
        self._load_ff_widgets ( )
        self._ff_loading = False
        self.ff_settings_toggle.connect ( 'toggled', lambda button: self._sync_side_panel ( ) )
        b ( 'ff_done_button' ).connect ( 'clicked', lambda button: self.ff_settings_toggle.set_active ( False ) )
        for radio in self.ff_radios.values ( ):
            radio.connect ( 'toggled', self.on_ff_setting_changed )
        for key, widget in self.ff_widgets.items ( ):
            signal = 'toggled' if key in ( "water", "bonded" ) else 'value-changed' if "cm5_scale" in key else 'changed'
            widget.connect ( signal, self.on_ff_setting_changed )

    def ff_settings ( self ):
        """ The Force Field Settings (defaults filled in). """
        settings = dict ( FF_SETTINGS_DEFAULTS )
        saved = self.vm_session.vm_config.gl_parameters.get ( 'builder_ff_settings' ) or { }
        settings.update ( { k: v for k, v in saved.items ( ) if k in FF_SETTINGS_DEFAULTS } )
        return settings

    def _load_ff_widgets ( self ):
        settings = self.ff_settings ( )
        self.ff_radios.get ( settings["force_field"], self.ff_radios["dyff"] ).set_active ( True )
        self.ff_widgets["method"].set_active_id ( settings["method"] )
        self.ff_widgets["water"].set_active ( bool ( settings["water"] ) )
        self.ff_widgets["bonded"].set_active ( bool ( settings["bonded"] ) )
        self.ff_widgets["cm5_scale"].set_value ( float ( settings["cm5_scale"] ) )
        self.ff_widgets["nb_model"].set_active_id ( settings["nb_model"] )
        self.ff_widgets["on_close"].set_active_id ( settings["on_close"] )
        self.ff_widgets["opls_cm5_scale"].set_value ( float ( settings["opls_cm5_scale"] ) )
        self._update_ff_panel_state ( settings )

    def _update_ff_panel_state ( self, settings ):
        self.builder.get_object ( 'builder_ff_description_label' ).set_text ( FF_DESCRIPTIONS[settings["force_field"]] )
        self.builder.get_object ( 'builder_ff_xtb_frame' ).set_sensitive ( settings["force_field"] == "dyff_xtb" )
        self.builder.get_object ( 'builder_ff_opls_frame' ).set_sensitive ( settings["force_field"] == "opls" )
        # . tblite (g-xTB) has no implicit solvent
        self.ff_widgets["water"].set_sensitive ( settings["method"] == "gfn2" )
        self.ff_widgets["nb_model"].set_sensitive ( settings["force_field"] != "none" )
        self.ff_widgets["on_close"].set_sensitive ( settings["force_field"] != "none" )

    def on_ff_setting_changed ( self, widget ):
        if self._ff_loading or ( isinstance ( widget, Gtk.RadioButton ) and not widget.get_active ( ) ):
            return
        settings = {
            "force_field": next ( key for key, radio in self.ff_radios.items ( ) if radio.get_active ( ) ),
            "method":      self.ff_widgets["method"].get_active_id ( ) or "gfn2",
            "water":       self.ff_widgets["water"].get_active ( ),
            "bonded":      self.ff_widgets["bonded"].get_active ( ),
            "cm5_scale":   round ( self.ff_widgets["cm5_scale"].get_value ( ), 2 ),
            "nb_model":    self.ff_widgets["nb_model"].get_active_id ( ) or "cutoff",
            "on_close":    self.ff_widgets["on_close"].get_active_id ( ) or "ask",
            "opls_cm5_scale": round ( self.ff_widgets["opls_cm5_scale"].get_value ( ), 2 ),
        }
        self._update_ff_panel_state ( settings )
        vm_config = self.vm_session.vm_config
        vm_config.gl_parameters['builder_ff_settings'] = settings
        try:
            vm_config.save_easyhybrid_config ( )
        except Exception as error:
            print ( "Force Field Settings not saved:", error )

    def _transform_selected_atom_ids ( self ):
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return [ ]
        selection = self.vm_session.selections.get ( self.vm_session.current_selection )
        selected_atoms = getattr ( selection, "selected_atoms", None ) or set ( )
        return sorted ( atom.atom_id for atom in selected_atoms
                        if atom.vm_object is target_object and target_object.atoms.get ( atom.atom_id ) is atom )

    def _refresh_transform_header ( self ):
        n = len ( self._transform_selected_atom_ids ( ) )
        text = "{} atom(s) selected".format ( n ) if n else "No atoms selected."
        if self.transform_header_label.get_text ( ) != text:
            self.transform_header_label.set_text ( text )

    def _enter_transform_mode ( self ):
        if self.placemode_toggle.get_active ( ):
            self.placemode_toggle.set_active ( False )
            self._transform_paused_editing = True
        selection_box = getattr ( self.vm_session, "selection_box_frame", None )
        if selection_box is not None:
            selection_box.change_toggle_button_selecting_mode_status ( False )
        self.vm_session.picking_selection_mode = False
        self._reset_transform_sliders ( )
        self._refresh_transform_header ( )
        self.vm_session.vm_glcore.queue_draw ( )
        if getattr ( self.main, "statusbar_main", None ) is not None:
            self.main.statusbar_main.push ( 1, "Transform: select atoms (click / shift+click / shift+drag), then use the sliders." )

    def _leave_transform_mode ( self ):
        if self._transform_paused_editing:
            self._transform_paused_editing = False
            self.placemode_toggle.set_active ( True )

    def _reset_transform_sliders ( self ):
        self._transform_last_values = { key: 0.0 for key in self.transform_scales }
        self._transform_gesture_atoms = None
        self._transform_suppress = True
        for scale in self.transform_scales.values ( ):
            scale.set_value ( 0.0 )
        self._transform_suppress = False

    def on_transform_slider_changed ( self, scale, key ):
        if self._transform_suppress:
            return
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return
        new_value = scale.get_value ( )
        delta = new_value - self._transform_last_values[key]
        self._transform_last_values[key] = new_value
        if abs ( delta ) < 1e-9:
            return
        from gui.windows.builder.atom_ops import ( translate_atoms, rotate_atoms_around_point,
                                                   selection_centroid, push_undo_snapshot )
        if self._transform_gesture_atoms is None:            # first move of this drag
            atom_ids = self._transform_selected_atom_ids ( )
            if not atom_ids:
                if getattr ( self.main, "statusbar_main", None ) is not None:
                    self.main.statusbar_main.push ( 1, "Transform: select some atoms first." )
                return
            push_undo_snapshot ( target_object )
            self._transform_gesture_atoms = atom_ids
        atom_ids = [ i for i in self._transform_gesture_atoms if i in target_object.atoms ]
        if key[0] == "t":
            d = [ 0.0, 0.0, 0.0 ]
            d["xyz".index ( key[1] )] = delta
            translate_atoms ( target_object, atom_ids, *d )
        else:
            import math
            axis = tuple ( 1.0 if a == key[1] else 0.0 for a in "xyz" )
            pivot = selection_centroid ( target_object, atom_ids )
            rotate_atoms_around_point ( target_object, atom_ids, pivot, axis, math.radians ( delta ) )
        # keyboard changes have no button release: sync once they stop
        from gi.repository import GLib
        if self._transform_sync_timer is not None:
            GLib.source_remove ( self._transform_sync_timer )
        self._transform_sync_timer = GLib.timeout_add ( 400, self._sync_after_transform, target_object )

    def _sync_after_transform ( self, target_object ):
        self._transform_sync_timer = None
        if target_object in self.vm_session.vm_objects_dic.values ( ):
            from gui.windows.builder.empty_object import sync_pdynamo_system
            sync_pdynamo_system ( target_object )
        return False

    def on_transform_slider_released ( self, widget, event ):
        """ End of a drag: the sliders go back to 0 so the next drag starts
        from a clean baseline (the pDynamo sync is the debounced one above). """
        from gi.repository import GLib
        GLib.idle_add ( lambda: ( self._reset_transform_sliders ( ), False )[1] )
        return False

    def _bind_delete_panel ( self ):
        """ [EN] 2026-10-17, user's request: what a Delete click removes --
        the atom, or its whole residue / molecule (bonded fragment) / chain
        (vm_session.builder_delete_scope, read by click_mode.handle_click_to_
        delete_atom()). """
        self.delete_scope_radios = { }
        for scope in ( "atom", "residue", "molecule", "chain" ):
            radio = self.builder.get_object ( 'delete_scope_{}_radio'.format ( scope ) )
            self.delete_scope_radios[scope] = radio
            radio.connect ( 'toggled', self.on_delete_scope_toggled, scope )

    def on_delete_scope_toggled ( self, button, scope ):
        if button.get_active ( ):
            self.vm_session.builder_delete_scope = scope

    def _bind_dihedral_panel ( self ):
        """ [EN] 2026-10-17, user's request: the dihedral controls are a page
        of the side panel (no Dihedral Angle window). click_mode.handle_click_
        to_pick_dihedral_atom() still calls arm/disarm_dihedral_slider(). """
        self.dihedral_adjustment   = self.builder.get_object ( 'builder_dihedral_adjustment' )
        self.dihedral_status_label = self.builder.get_object ( 'builder_dihedral_status_label' )
        self.dihedral_controls     = [ self.builder.get_object ( name ) for name in
                                       ( 'builder_dihedral_scale', 'builder_dihedral_spin', 'builder_dihedral_preset_grid' ) ]
        self._suppress_dihedral_signal = False
        self._dihedral_sync_timer = None
        self.dihedral_adjustment.connect ( 'value-changed', self.on_dihedral_value_changed )
        for name, value in ( ( 'dihedral_preset_180', 180.0 ), ( 'dihedral_preset_p60', 60.0 ),
                             ( 'dihedral_preset_m60', -60.0 ), ( 'dihedral_preset_0', 0.0 ) ):
            self.builder.get_object ( name ).connect ( 'clicked', lambda button, v: self.dihedral_adjustment.set_value ( v ), value )
        self._set_dihedral_controls_sensitive ( False )

    def _set_dihedral_controls_sensitive ( self, sensitive ):
        for widget in self.dihedral_controls:
            if widget is not None:
                widget.set_sensitive ( sensitive )

    def on_dihedral_value_changed ( self, adjustment ):
        """ Rotates the armed dihedral (vm_session.builder_dihedral_axis/
        _subgroup) so that its angle -- re-measured every call, never
        accumulated -- becomes the new value (moved unchanged from the old
        DihedralAngleWindow). The pDynamo system is synced once the value
        stops changing (debounced), not on every slider motion. """
        if self._suppress_dihedral_signal:
            return
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        axis          = getattr ( self.vm_session, "builder_dihedral_axis", None )
        subgroup      = getattr ( self.vm_session, "builder_dihedral_subgroup", None )
        picked        = getattr ( self.vm_session, "builder_dihedral_pick_atoms", None )
        if target_object is None or axis is None or subgroup is None or not picked or len ( picked ) < 4:
            return

        import math
        from util.geometric_analysis import get_dihedral
        try:
            current_angle = get_dihedral ( target_object, *[ a.atom_id for a in picked ] )
        except ValueError:
            return
        delta_deg = adjustment.get_value ( ) - current_angle
        while delta_deg > 180.0:
            delta_deg -= 360.0
        while delta_deg <= -180.0:
            delta_deg += 360.0
        if abs ( delta_deg ) < 1e-9:
            return

        from gui.windows.builder.atom_ops import rotate_atoms_around_bond
        rotate_atoms_around_bond ( target_object, axis[0], axis[1], subgroup, math.radians ( delta_deg ) )
        self._set_dihedral_status ( adjustment.get_value ( ) )

        from gi.repository import GLib
        if self._dihedral_sync_timer is not None:
            GLib.source_remove ( self._dihedral_sync_timer )
        self._dihedral_sync_timer = GLib.timeout_add ( 400, self._sync_after_dihedral_change, target_object )

    def _sync_after_dihedral_change ( self, target_object ):
        self._dihedral_sync_timer = None
        if target_object in self.vm_session.vm_objects_dic.values ( ):
            from gui.windows.builder.empty_object import sync_pdynamo_system
            sync_pdynamo_system ( target_object )
        return False

    def _start_session ( self ):

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
        self.ff_settings_toggle.set_active ( False )
        self.placemode_toggle.set_label ( "Editing: ON" )
        self.tool_add_radio.set_active ( True )
        self.add_notebook.set_current_page ( 0 )      # "Chemical"
        self.element_C_radio.set_active ( True )
        self.bond_order_single_radio.set_active ( True )
        self.clean_up_adjust_hydrogens_checkbutton.set_active ( True )
        self.vm_session.builder_delete_scope = "atom"
        self.delete_scope_radios["atom"].set_active ( True )
        self._chemical_state = None
        self._selected_structure = None
        self.vm_session.builder_selected_structure = None
        self.structure_status_label.set_text ( "No structure selected." )
        if self.structure_preview is not None:
            self.structure_preview.clear ( )
        self._populate_structure_tree ( )

        self.toolbar_edit.set_no_show_all ( False )
        self.toolbar_edit.show_all ( )
        self.visible = True
        self._sync_side_panel ( )
        self.sync_main_toggle ( )
        # keep Undo/Redo greyed out when there is nothing to undo/redo --
        # edits happen in click_mode.py/vismol_glcore.py, so poll cheaply
        from gi.repository import GLib
        self._undo_redo_timer = GLib.timeout_add ( 300, self._update_undo_redo_sensitivity )
        self._update_undo_redo_sensitivity ( )

    def sync_main_toggle ( self ):
        """ Main toolbar's "Edit and Build" toggle (_show_cell) follows the
        real session state -- the session can also be opened from the
        treeview ("Edit in Builder") or closed from toolbar_edit's Close. """
        button = self.main.builder.get_object ( '_show_cell' )
        if button is None or button.get_active ( ) == self.visible:
            return
        self.main._syncing_builder_toggle = True
        try:
            button.set_active ( self.visible )
        finally:
            self.main._syncing_builder_toggle = False

    def _sync_side_panel ( self ):
        """ [EN] 2026-10-17: the right side panel shows what the active tool
        needs -- elements + fragments with "Add", the bond order radios with
        "Bond Order" (used to be a row shown/hidden inside the old window);
        no panel for the other tools. """
        if not self.visible:
            return
        if self.ff_settings_toggle.get_active ( ):
            self.main.show_sidebar_page ( "builder_ff_settings", "Force Field Settings" )
        elif self.tool_add_radio.get_active ( ):
            self.main.show_sidebar_page ( "builder_add", "Add" )
        elif self.tool_bond_order_radio.get_active ( ):
            self.main.show_sidebar_page ( "builder_bond_order", "Bond order" )
        elif self.tool_delete_radio.get_active ( ):
            self.main.show_sidebar_page ( "builder_delete", "Delete" )
        elif self.tool_rotate_dihedral_radio.get_active ( ):
            self.main.show_sidebar_page ( "builder_dihedral", "Dihedral" )
        elif self.transform_selection_button.get_active ( ):
            self.main.show_sidebar_page ( "builder_transform", "Transform selection" )
        else:
            self.main.hide_sidebar ( )

    def _update_undo_redo_sensitivity ( self ):
        if not self.visible:
            self._undo_redo_timer = None
            return False
        if self.transform_selection_button.get_active ( ):
            self._refresh_transform_header ( )
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        can_undo = bool ( getattr ( target_object, "undo_stack", None ) )
        can_redo = bool ( getattr ( target_object, "redo_stack", None ) )
        if self.undo_button is not None and self.undo_button.get_sensitive ( ) != can_undo:
            self.undo_button.set_sensitive ( can_undo )
        if self.redo_button is not None and self.redo_button.get_sensitive ( ) != can_redo:
            self.redo_button.set_sensitive ( can_redo )
        return True

    def close_window ( self, *args, allow_cancel = False ):
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
        # [EN] 2026-10-04 BUG FIX: the object may have been deleted while the
        # Builder was open (treeview Delete -> main_window._on_vobjects_removed());
        # then there is nothing to ask, keep, discard or fold back -- answering
        # "Yes" used to resurrect it as an orphan pDynamo system.
        if target_object is not None and target_object not in self.vm_session.vm_objects_dic.values ( ):
            target_object = None

        # [EN] 2026-10-04, user's request: when the session leaves a system
        # behind, ask whether to incorporate it with DYFF and the Builder's
        # customizations (edited parameters, partial charges, atom types).
        # "No" keeps today's behaviour (a system without force field);
        # "Cancel" (only from the Close button / window X) keeps editing.
        incorporate = False
        if self._session_leaves_a_system ( target_object ):
            answer = self._ask_incorporate_with_dyff ( target_object, allow_cancel )
            if answer is None:
                return True
            incorporate = answer
        final_vobject = None

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
            final_vobject = finish_editing_existing_system ( target_object, carry_customizations = incorporate )
        elif target_object is not None:
            final_vobject = target_object
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
        # the Solvate window works on this session's object too
        solvate_window = getattr ( self.main, "solvate_window", None )
        if solvate_window is not None and solvate_window.visible:
            solvate_window.close_window ( )
        # the DYFF Parameters window (and its torsion fit window) is a Builder
        # tool on this session's object: it used to stay open on it
        dyff_window = getattr ( self.main, "dyff_parameters_window", None )
        if dyff_window is not None and dyff_window.visible:
            dyff_window.close_window ( )
        # [EN] The Atom Types window's own 3D-view index-label overlay
        # (see atom_types_window.py's _refresh_rows()) points at atoms
        # belonging to the object THIS Builder session was editing --
        # stale/dangling once that session ends, regardless of whether
        # the Atom Types window itself is still open or already closed.
        self.vm_session.builder_atom_types_labeled_atoms = [ ]
        timer = getattr ( self, "_undo_redo_timer", None )
        if timer is not None:
            from gi.repository import GLib
            GLib.source_remove ( timer )
            self._undo_redo_timer = None
        self.toolbar_edit.hide ( )
        self.toolbar_edit.set_no_show_all ( True )
        self.main.hide_sidebar ( )
        self.visible = False
        self.sync_main_toggle ( )
        if incorporate and final_vobject is not None:
            self._incorporate_with_dyff ( final_vobject, rebuild = final_vobject is target_object )
        return True

    def _session_leaves_a_system ( self, target_object ):
        """ True when closing keeps a system with atoms: a blank session that
            got real edits (promoted) or an "Edit in Builder" session. """
        if target_object is None or not getattr ( target_object, "atoms", None ):
            return False
        if getattr ( target_object, "builder_edit_source_e_id", None ) is not None:
            return True
        return not getattr ( target_object, "is_builder_only", False ) and getattr ( target_object, "e_id", None ) is not None

    def _ask_incorporate_with_dyff ( self, target_object, allow_cancel ):
        """ True (incorporate), False (plain system, as before) or None (cancel).
        [EN] 2026-10-17: follows the Force Field Settings -- "None" never
        asks, "Apply automatically" doesn't ask, the question names the
        chosen force field. """
        from gi.repository import Gtk
        from gui.windows.builder.dyff_parameters import builder_customizations, customizations_text
        library = getattr ( target_object, "dyff_library_options", None )
        settings = self.ff_settings ( )
        force_field = settings["force_field"]
        if not library and force_field == "none":
            return False
        if settings["on_close"] == "apply":
            return True
        dialog = Gtk.MessageDialog ( transient_for = getattr ( self, "window", None ) or self.main.window, modal = True,
                                     message_type = Gtk.MessageType.QUESTION, buttons = Gtk.ButtonsType.NONE,
                                     text = "Incorporate '{}' with the {} force field?".format (
                                         target_object.name, "DYFF" if library else FF_NAMES[force_field] ) )
        if library:
            # . [EN] 2026-10-08: a system prepared with DYFF + residue library (Prepare DYFF
            #   System): the same parametrization is applied again to the edited structure
            dialog.format_secondary_text (
                "This system was parametrized with DYFF + residue library. Yes: the library is applied again to the "
                "edited structure -- residues that still match it get the same parameters; the edited region (e.g. a "
                "metal site with a bonded hydroxide) is outside the library and is refined with {} (CM5 charges, "
                "Hessian bonds/angles, torsion library); metal types and charges are read from the structure.\n\n"
                "No: the system is kept without a force field.".format (
                    { "gxtb": "g-xTB" }.get ( library.get ( "method" ), "GFN2-xTB" ) ) )
        elif force_field == "dyff":
            lines = customizations_text ( builder_customizations ( target_object ) )
            dialog.format_secondary_text (
                "Yes: the system gets the DYFF MM model with the customizations made in the Builder:\n"
                + "".join ( "    \u2022 {}\n".format ( line ) for line in lines )
                + "\nNo: the system is kept without a force field (as before)."
                + "\n\n(Builder > Force Field Settings chooses the force field.)" )
        else:
            dialog.format_secondary_text (
                "Yes: {}\n\nNo: the system is kept without a force field.\n\n"
                "(Builder > Force Field Settings chooses the force field.)".format ( self._ff_settings_text ( settings ) ) )
        if allow_cancel:
            dialog.add_button ( "Cancel", Gtk.ResponseType.CANCEL )
        dialog.add_button ( "No, plain system", Gtk.ResponseType.NO )
        dialog.add_button ( "Yes, re-apply DYFF + residue library" if library else
                            "Yes, incorporate with {}".format ( FF_NAMES[force_field] ), Gtk.ResponseType.YES )
        dialog.set_default_response ( Gtk.ResponseType.YES )
        response = dialog.run ( )
        dialog.destroy ( )
        if response == Gtk.ResponseType.YES: return True
        if response == Gtk.ResponseType.NO:  return False
        return None if allow_cancel else False

    def _incorporate_with_dyff ( self, vismol_object, rebuild = True ):
        """ DYFF + customizations on the system the session left. rebuild: the
            pDynamo system is rebuilt from the object first (current
            coordinates; a fold-back has just done it). """
        from gui.windows.builder.dyff_parameters import incorporate_with_dyff
        from gui.windows.builder.empty_object import sync_pdynamo_system
        if rebuild:
            sync_pdynamo_system ( vismol_object )
        if getattr ( vismol_object, "dyff_library_options", None ):
            self.main.reapply_dyff_library ( vismol_object, extra_message = " (re-applied after editing in the Builder)" )
            return
        system = self.main.p_session.psystem.get ( getattr ( vismol_object, "e_id", None ) )
        # [EN] 2026-10-17: the force field chosen in Force Field Settings
        settings = self.ff_settings ( )
        force_field = settings["force_field"]
        if force_field == "none":
            return
        if force_field == "dyff_xtb":
            self._incorporate_with_dyff_xtb ( vismol_object, system, settings )
            return
        if force_field == "opls":
            self._incorporate_with_opls ( vismol_object, system, settings )
            return
        ok, message = incorporate_with_dyff ( self.main, vismol_object )
        if ok:
            message += self._apply_ff_nb_model ( system, settings )
        self.main.refresh_main_statusbar ( )
        if ok:
            self.main.bottom_notebook.status_teeview_add_new_item ( message = message, system = system )
        else:
            self.main.simple_dialog.error ( msg = message + "\nThe system was kept without a force field." )

    def _incorporate_with_dyff_xtb ( self, vismol_object, system, settings ):
        """ DYFF + xTB refinement: the residue library pipeline of Prepare
            DYFF System (main_window._run_dyff_library) -- for a drawn
            molecule no residue matches, so every atom gets CM5 charges,
            Hessian bonds/angles and torsions. Runs in a worker thread; the
            options stay on the object (dyff_library_options), so a later
            Edit in Builder re-applies them. """
        from pdynamo import residue_library as RL
        if system is None:
            self.main.simple_dialog.error ( msg = "No pDynamo system for '{}'.".format ( vismol_object.name ) )
            return
        if getattr ( system, "qcModel", None ) is not None:
            self.main.simple_dialog.error ( msg = "'{}' has a QC model: DYFF was not assigned.".format ( system.label ) )
            return
        try:
            plan = RL.plan ( system )
        except Exception as error:
            self.main.simple_dialog.error ( msg = "Residue library: {}\nThe system was kept without a force field.".format ( error ) )
            return
        self.main._run_dyff_library ( system, vismol_object, plan, use_water = settings["water"] and settings["method"] == "gfn2",
                                      use_bonded = settings["bonded"], method = settings["method"],
                                      cm5_scale = float ( settings["cm5_scale"] ),
                                      extra_message = " (Builder Force Field Settings)",
                                      on_success = lambda: self._apply_ff_nb_model ( system, settings ) )

    def _incorporate_with_opls ( self, vismol_object, system, settings ):
        """ [EN] 2026-10-17, user's request: the OPLS-AA set as it is; only when
        parameters are missing, a dialog offers the Prepare OPLS System procedure to
        create them (no solvation -- a separate, later step). See
        gui/windows/setup/opls_assign.py. A Builder molecule keeps its atoms
        through that procedure, so the object stays the same. """
        from gui.windows.setup.opls_assign import assign_opls
        if system is None:
            self.main.simple_dialog.error ( msg = "No pDynamo system for '{}'.".format ( vismol_object.name ) )
            return
        assign_opls ( self.main, system, cm5_scale = float ( settings["opls_cm5_scale"] ),
                      after = lambda final_system: self._apply_ff_nb_model ( final_system, settings ) )

    def _apply_ff_nb_model ( self, system, settings ):
        """ define_MMModel() already gives a cut-off NB model when there is
            none; "Full" replaces it. Returns a note for the message. """
        if system is None or settings["nb_model"] != "full":
            return ""
        if self.main.p_session.define_NBModel ( _type = 0, system = system ):
            return " Nonbonded model: full (no cut-off)."
        return " (The full nonbonded model could not be defined; cut-off kept.)"

    def _ff_settings_text ( self, settings ):
        force_field = settings["force_field"]
        if force_field == "dyff_xtb":
            text = ( "DYFF with CM5 charges (scale {:.2f}){} computed with {}{}; runs in the background.".format (
                     float ( settings["cm5_scale"] ),
                     ", Hessian bonds/angles and torsions" if settings["bonded"] else "",
                     { "gxtb": "g-xTB" }.get ( settings["method"], "GFN2-xTB" ),
                     " in implicit water" if settings["water"] and settings["method"] == "gfn2" else "" ) )
        elif force_field == "opls":
            text = ( "OPLS-AA. If parameters are missing, you will be asked whether to create them with the "
                     "Prepare OPLS System procedure (CM5 scale {:.2f}).".format ( float ( settings["opls_cm5_scale"] ) ) )
        else:
            text = "DYFF ({}).".format ( FF_DESCRIPTIONS["dyff"] )
        if settings["nb_model"] == "full":
            text += " Nonbonded model: full (no cut-off)."
        return text

    # ------------------------------------------------------------------
    #  Signal handlers (referenced by name in builder_sidebar.glade)
    # ------------------------------------------------------------------

    def on_placemode_toggled ( self, button ):
        """ The lighter-weight PAUSE toggle -- see this module's design
        note at the top for how this differs from actually closing the
        window. """
        is_active = button.get_active ( )
        self.vm_session.builder_atom_mode = is_active
        if is_active:
            self._transform_paused_editing = False   # turned back on by hand
        button.set_label ( "Editing: ON" if is_active else "Editing: OFF" )
        if self.visible and getattr ( self.main, "statusbar_main", None ) is not None:
            self.main.statusbar_main.push ( 1, "Builder editing ON." if is_active else
                                            "Builder editing paused (OFF): clicks only rotate/inspect the view." )

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
            self.vm_session.builder_tool = "add_structure" if self._structures_tab_active ( ) else "add"
        elif button is self.tool_delete_radio:
            self.vm_session.builder_tool = "delete"
        elif button is self.tool_move_radio:
            self.vm_session.builder_tool = "move"
        elif button is self.tool_bond_order_radio:
            self.vm_session.builder_tool = "bond_order"
        elif button is self.tool_rotate_dihedral_radio:
            self.vm_session.builder_tool = "rotate_dihedral"
        elif button is self.transform_selection_button:
            self.vm_session.builder_tool = "transform"
            self._enter_transform_mode ( )
        if self.vm_session.builder_tool != "transform":
            self._leave_transform_mode ( )

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
        if self.ff_settings_toggle.get_active ( ):
            self.ff_settings_toggle.set_active ( False )   # its handler syncs the panel
        else:
            self._sync_side_panel ( )

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

        # same toggle delivered twice (glade signal + the Python connection
        # made in open_window()) -> nothing new to do
        last_button, last_fragment = self._last_quick_fragment
        if ( button is last_button and last_fragment is not None
                and self.vm_session.builder_selected_fragment is last_fragment
                and getattr ( self.vm_session, "builder_tool", None ) == "attach_fragment" ):
            return

        file_path = os.path.join ( self.main.home, "src/gui/windows/builder/fragments", filename )
        try:
            fragment = load_fragment ( file_path )
        except FragmentError as error:
            self.main.simple_dialog.error ( msg = str ( error ) )
            return
        self._last_quick_fragment = ( button, fragment )

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
        on_tool_changed()/on_element_changed() above. [EN] 2026-10-03:
        this row is shown only with the "Bond Order" tool and no longer
        affects "Add" (a dragged bond's order comes from the drag
        distance). It still applies to the 'b' key, to
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

        # [EN] 2026-10-04: opens the DYFF Parameters window (replaces Atom Types;
        # see dyff_parameters_window.py) -- on the selected atoms, or on the whole
        # molecule when nothing is selected.
        if not target_object.atoms:
            if getattr ( self.main, "statusbar_main", None ) is not None:
                self.main.statusbar_main.push ( 1, "DYFF Parameters: the molecule has no atoms yet." )
            return
        self.main.dyff_parameters_window.open_window ( target_object, atom_ids )

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
        radio), so no change needed at that call site.
        [EN] 2026-10-17: drives the Dihedral page of the side panel instead
        (user's request: no auxiliary window). """
        self._suppress_dihedral_signal = True
        self.dihedral_adjustment.set_value ( angle_deg )
        self._suppress_dihedral_signal = False
        self._set_dihedral_controls_sensitive ( True )
        self._set_dihedral_status ( angle_deg )

    def _set_dihedral_status ( self, angle_deg ):
        picked = getattr ( self.vm_session, "builder_dihedral_pick_atoms", None ) or [ ]
        names = "-".join ( a.name for a in picked[:4] )
        self.dihedral_status_label.set_text ( "{}: {:.1f}\u00b0".format ( names or "Dihedral", angle_deg ) )

    def disarm_dihedral_slider ( self ):
        """ Counterpart to arm_dihedral_slider() above (tool switched away,
        pick cancelled, or the Builder closing): controls disabled until
        the next 4-atom pick. """
        if not self._widgets_bound:
            return
        self._set_dihedral_controls_sensitive ( False )
        self.dihedral_status_label.set_text ( "No dihedral picked." )

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
        self._update_undo_redo_sensitivity ( )

    def on_redo_button_clicked ( self, button ):
        """ [EN] 2026-10-03: re-applies the last action undone by Undo
        (atom_ops.redo()); any new edit clears what can be redone. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return
        atom_ops_redo ( target_object )
        self._update_undo_redo_sensitivity ( )

    def _create_add_hydrogens_button ( self ):
        existing = self.builder.get_object ( "add_hydrogens_button" )
        if existing is not None:          # already moved into the glade
            return existing
        anchor = self.optimize_dyff_button
        parent = anchor.get_parent ( ) if anchor is not None else None
        if not isinstance ( parent, Gtk.Box ):
            return None
        button = Gtk.Button ( label = "Add Hydrogens..." )
        button.set_tooltip_text ( "Adds the missing hydrogens of the molecule being edited. Standard residues "
                                  "use pDynamo's residue templates (names HA, HB2, ...; protonation reviewed in a "
                                  "dialog), other molecules are protonated by OpenBabel. Every new hydrogen is "
                                  "placed inside its own residue. Can be undone." )
        button.connect ( "clicked", self.on_add_hydrogens_button_clicked )
        children = parent.get_children ( )
        parent.pack_start ( button, False, False, 0 )
        if anchor in children:
            parent.reorder_child ( button, children.index ( anchor ) + 1 )
        # copy the anchor's margins so it lines up with its neighbours
        for side in ( "start", "end", "top", "bottom" ):
            getattr ( button, "set_margin_" + side ) ( getattr ( anchor, "get_margin_" + side ) ( ) )
        button.show ( )
        return button

    def on_add_hydrogens_button_clicked ( self, button ):
        """ [EN] 2026-10-02: hydrogens for the molecule being edited, through
        the SAME residue-template machinery as the main "Add Missing
        Hydrogens" tool (util/hydrogen_builder.py: pDynamo PDB component
        templates for standard residues, OpenBabel for the rest), opened on
        this object with its protonation-review dialog. The Builder's linked
        pDynamo system now carries chains/residues (empty_object.
        _build_sequenced_pdynamo_system()), which is what lets the templates
        recognise each residue; the rebuilt system keeps every hydrogen
        INSIDE its residue. Before each rebuild the Builder's own undo
        snapshot is taken; afterwards bond orders/charges are re-derived
        (atom_ops.rederive_bonding_after_rebuild()). """
        from gui.windows.builder.atom_ops import push_undo_snapshot, rederive_bonding_after_rebuild
        from gui.windows.builder.empty_object import sync_pdynamo_system
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None or not target_object.atoms:
            self.main.simple_dialog.info ( msg = "There is no molecule in the Builder yet." )
            return
        sync_pdynamo_system ( target_object )
        if getattr ( target_object, "e_id", None ) is None or self.main.p_session.psystem.get ( target_object.e_id ) is None:
            self.main.simple_dialog.error ( msg = "Could not build the pDynamo system for this molecule." )
            return

        def before ( vm_object ):
            push_undo_snapshot ( vm_object )

        def after ( vm_object ):
            # formal charges set by the rebuild (template variants, OpenBabel
            # protonation of ligands) live on the new pDynamo atoms -- copy
            # them onto the Builder's atoms before re-deriving bond orders
            system = self.main.p_session.psystem.get ( vm_object.e_id )
            if system is not None and len ( system.atoms ) == len ( vm_object.atoms ):
                for k, p_atom in enumerate ( system.atoms ):
                    vm_object.atoms[k].formal_charge = int ( getattr ( p_atom, "formalCharge", 0 ) or 0 )
            rederive_bonding_after_rebuild ( vm_object )
            sync_pdynamo_system ( vm_object )

        self.main.prepare_add_hydrogens_window.open_window (
                target_vobject = target_object, on_before_apply = before, on_after_apply = after )

    def on_solvate_button_clicked ( self, button ):
        """ [EN] 2026-10-03: opens the Solvate window (solvate_window.py) on
        the molecule being edited. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None or not target_object.atoms:
            self.main.simple_dialog.info ( msg = "There is no molecule in the Builder yet." )
            return
        self.main.solvate_window.open_window ( target_object )

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
        return self.close_window ( allow_cancel = True )
