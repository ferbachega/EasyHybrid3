#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: xtb_fragment_charges_window
#
#  Copyright 2022-2025 Fernando Bachega
#
"""
xtb_fragment_charges_window
===========================

GTK window (Phase 2) for the self-consistent fragment-based MM charge tool
(engine in util/xtb_fragment_charges.py).

Usage flow:
  1. choose the level (residue/segment/chain) -> "Build fragments" populates the
     treeview with SUGGESTED formal charge and multiplicity;
  2. the user can edit charge/multiplicity of each fragment in the treeview
     (the treeview is the source of truth before running);
  3. choose method (GFN0/1/2), charge model (CM5/Mulliken), tolerance,
     boundary mode, xTB path and number of processes;
  4. "Run" -> runs the self-consistent loop in a separate THREAD (each fragment
     is a process, via multiprocessing.Pool); the progress bar updates per cycle;
  5. "Apply to system" -> writes the final charges into system.mmState.charges.

Integration: instantiate with (main, system) and call open_window(). E.g.:
    win = XtbFragmentChargesWindow(main=self.main, system=active_system)
    win.open_window()
"""

import os
import threading

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

from util import xtb_fragment_charges as xfc


# --------------------------------------------------------------------------- #
#  Clustal X color palette (by amino-acid type), RGB 0-255.                    #
#  Simplified "per-residue" version (without the column-dependent alignment    #
#  rules): each residue group gets a fixed color.                              #
#  [VERIFY: confirm these are the Clustal colors you want for the paper]        #
# --------------------------------------------------------------------------- #
CLUSTAL_COLORS = {
    # hydrophobic / blue
    "A": (128, 160, 240), "I": (128, 160, 240), "L": (128, 160, 240),
    "M": (128, 160, 240), "F": (128, 160, 240), "W": (128, 160, 240),
    "V": (128, 160, 240), "C": (128, 160, 240),
    # red
    "K": (240, 21, 5), "R": (240, 21, 5),
    # magenta
    "D": (192, 72, 192), "E": (192, 72, 192),
    # green
    "N": (21, 194, 194), "Q": (21, 194, 194),
    "S": (21, 194, 194), "T": (21, 194, 194),
    # orange
    "G": (240, 144, 72),
    # yellow
    "P": (192, 192, 0),
    # cyan
    "H": (21, 164, 164), "Y": (21, 164, 164),
}
DEFAULT_RESIDUE_COLOR = (200, 200, 200)  # gray for anything that does not match

# 3-letter -> 1-letter code (to match the palette from the residue name)
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
    # common protonation variants -> same base residue
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "CYM": "C",
    "ASH": "D", "GLH": "E", "LYN": "K",
}


def _residue_one_letter(fragment_key):
    """Extract the 1-letter residue code from the 'CHAIN/RESNAME/SEQ' key."""
    try:
        resname = fragment_key.split("/")[1].upper()
    except Exception:
        return None
    if len(resname) == 1:
        return resname
    return _THREE_TO_ONE.get(resname[:3])


def clustal_rgb(fragment_key):
    """Clustal RGB color (0-255) for the fragment key; gray if unknown."""
    one = _residue_one_letter(fragment_key)
    return CLUSTAL_COLORS.get(one, DEFAULT_RESIDUE_COLOR)


def make_color_swatch(rgb, width=18, height=18):
    """Create a GdkPixbuf.Pixbuf filled with the rgb color (0-255)."""
    from gi.repository import GdkPixbuf
    r, g, b = int(rgb[0]), int(rgb[1]), int(rgb[2])
    a = 255
    swatch = GdkPixbuf.Pixbuf.new(GdkPixbuf.Colorspace.RGB, True, 8, width, height)
    swatch.fill((r << 24) | (g << 16) | (b << 8) | a)  # RGBA
    return swatch


# colunas do liststore:
# color(pixbuf), include(bool), key, natoms, formal_charge, multiplicity, status
COL_COLOR, COL_INCLUDE, COL_KEY, COL_NAT, COL_CHG, COL_MULT, COL_STATUS = range(7)


class XtbFragmentChargesWindow:
    """Window for the fragment-based charge tool (xTB)."""

    def __init__(self, main=None, system=None):
        self.main = main.window
        self.easy_session = main
        self.system = system
        self.window = None
        self.fragments = []          # lista de xfc.Fragment
        self.result = None           # dict returned by the engine
        self._running = False
        self._selection_counter = 0  # numera os fragmentos importados de selecao

    # ------------------------------------------------------------------ #
    #  Construcao da UI                                                    #
    # ------------------------------------------------------------------ #
    def _on_choose_tmpdir(self, _button):
        """Open a FOLDER chooser for the xTB temporary files."""
        dialog = Gtk.FileChooserDialog(
            title="Select temporary folder",
            transient_for=self.window,
            action=Gtk.FileChooserAction.SELECT_FOLDER)
        dialog.add_buttons(
            "Cancel", Gtk.ResponseType.CANCEL,
            "Select", Gtk.ResponseType.OK)
        if dialog.run() == Gtk.ResponseType.OK:
            self.entry_tmpdir.set_text(dialog.get_filename())
        dialog.destroy()

    def _system_color_swatch(self):
        """Pixbuf of the small color square that identifies the system.

        Reuses the SAME function as the main window
        (get_colorful_square_pixel_buffer), so the square matches exactly the one
        shown in the main treeview. The color comes from system.e_color_palette['C'].
        Returns None if not possible (import failed or system without palette),
        in which case the window shows text only.
        """
        if self.system is None:
            return None
        try:
            from gui.windows.setup.setup_interface import get_colorful_square_pixel_buffer
            return get_colorful_square_pixel_buffer(self.system)
        except Exception:
            # fallback: build directly from the palette, without relying on the import
            try:
                from gi.repository import GdkPixbuf
                color = self.system.e_color_palette["C"]
                r, g, b = (int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
                sw = GdkPixbuf.Pixbuf.new(GdkPixbuf.Colorspace.RGB, True, 8, 20, 20)
                sw.fill((r << 24) | (g << 16) | (b << 8) | 255)
                return sw
            except Exception:
                return None

    def _system_name(self):
        """Name/identification of the active system, shown at the top of the window.

        Uses system.label (the system name in EasyHybrid); adds e_tag when
        available. Falls back to a generic text if nothing exists.
        """
        if self.system is None:
            return "(no system)"
        label = getattr(self.system, "label", None)
        tag = getattr(self.system, "e_tag", None)
        if label and tag:
            return "{} ({})".format(label, tag)
        if label:
            return str(label)
        return "(unnamed system)"

    def open_window(self):
        if self.window is not None:
            self.window.present()
            return

        self.window = Gtk.Window(title="xTB fragment charges")
        self.window.set_default_size(640, 560)
        self.window.set_border_width(8)
        self.window.connect("destroy", self._on_destroy)

        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
        self.window.add(vbox)

        # ---- active system identification (color swatch + name) ----
        # shows which system is being worked on, with the SAME small color
        # square used in the main-window treeview, for quick visual
        # recognition when several systems are loaded.
        row_sys = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        swatch = self._system_color_swatch()
        if swatch is not None:
            self.img_system = Gtk.Image.new_from_pixbuf(swatch)
            row_sys.pack_start(self.img_system, False, False, 0)
        self.label_system = Gtk.Label()
        self.label_system.set_xalign(0.0)
        self.label_system.set_markup(
            "<b>System:</b> {}".format(self._system_name()))
        row_sys.pack_start(self.label_system, False, False, 0)
        vbox.pack_start(row_sys, False, False, 0)

        # ---- row 1: level + build / import selection ----
        row1 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        row1.pack_start(Gtk.Label(label="Fragment level:"), False, False, 0)
        self.combo_level = Gtk.ComboBoxText()
        for lv in ("residue", "segment", "chain", "selection"):
            self.combo_level.append_text(lv)
        self.combo_level.set_active(0)
        self.combo_level.connect("changed", self._on_level_changed)
        row1.pack_start(self.combo_level, False, False, 0)

        self.btn_build = Gtk.Button(label="Build fragments")
        self.btn_build.connect("clicked", self.on_build_fragments)
        row1.pack_start(self.btn_build, False, False, 0)

        # button to import the current vismol selection as a fragment
        # (visible only in 'selection' mode)
        self.btn_import_sel = Gtk.Button(label="Import selection fragment")
        self.btn_import_sel.connect("clicked", self.on_import_selection)
        row1.pack_start(self.btn_import_sel, False, False, 0)
        vbox.pack_start(row1, False, False, 0)

        # ---- fragment treeview (editable) ----
        # store: color(pixbuf), include(bool), key, natoms, charge, mult, status
        from gi.repository import GdkPixbuf
        self.store = Gtk.ListStore(GdkPixbuf.Pixbuf, bool, str, int, int, int, str)
        self.treeview = Gtk.TreeView(model=self.store)

        # column 0: color swatch (Clustal)
        renderer_color = Gtk.CellRendererPixbuf()
        col_color = Gtk.TreeViewColumn("", renderer_color, pixbuf=COL_COLOR)
        self.treeview.append_column(col_color)

        # coluna 1: checkbox incluir/ignorar
        renderer_toggle = Gtk.CellRendererToggle()
        renderer_toggle.connect("toggled", self._on_include_toggled)
        col_toggle = Gtk.TreeViewColumn("Use", renderer_toggle, active=COL_INCLUDE)
        self.treeview.append_column(col_toggle)

        self._add_text_column("Fragment", COL_KEY, editable=False)
        self._add_text_column("Atoms", COL_NAT, editable=False)
        self._add_text_column("Charge", COL_CHG, editable=True)
        self._add_text_column("Mult", COL_MULT, editable=True)
        self._add_text_column("Status", COL_STATUS, editable=False)

        scroll = Gtk.ScrolledWindow()
        scroll.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        scroll.add(self.treeview)
        vbox.pack_start(scroll, True, True, 0)

        # context menu (right-click) to (de)select all
        self.treeview.connect("button-press-event", self._on_treeview_button_press)

        # ---- calculation options ----
        grid = Gtk.Grid(column_spacing=8, row_spacing=6)

        grid.attach(Gtk.Label(label="Method:"), 0, 0, 1, 1)
        self.combo_method = Gtk.ComboBoxText()
        for m in ("gfn1", "gfn2", "gfn0"):
            self.combo_method.append_text(m)
        self.combo_method.set_active(0)
        grid.attach(self.combo_method, 1, 0, 1, 1)

        grid.attach(Gtk.Label(label="Charge model:"), 2, 0, 1, 1)
        self.combo_charge = Gtk.ComboBoxText()
        for c in ("cm5", "mulliken"):
            self.combo_charge.append_text(c)
        self.combo_charge.set_active(0)
        grid.attach(self.combo_charge, 3, 0, 1, 1)

        grid.attach(Gtk.Label(label="Boundary:"), 0, 1, 1, 1)
        self.combo_boundary = Gtk.ComboBoxText()
        for b in ("redistribute", "discard", "keep_on_boundary"):
            self.combo_boundary.append_text(b)
        self.combo_boundary.set_active(0)
        grid.attach(self.combo_boundary, 1, 1, 1, 1)

        grid.attach(Gtk.Label(label="Tolerance:"), 2, 1, 1, 1)
        self.entry_tol = Gtk.Entry()
        self.entry_tol.set_text("0.01")
        self.entry_tol.set_width_chars(8)
        grid.attach(self.entry_tol, 3, 1, 1, 1)

        grid.attach(Gtk.Label(label="Processes:"), 0, 2, 1, 1)
        self.entry_nproc = Gtk.Entry()
        self.entry_nproc.set_text("4")
        self.entry_nproc.set_width_chars(8)
        grid.attach(self.entry_nproc, 1, 2, 1, 1)
        

        grid.attach(Gtk.Label(label="Max cycles:"), 2, 2, 1, 1)
        self.entry_maxcyc = Gtk.Entry()
        self.entry_maxcyc.set_text("25")
        self.entry_maxcyc.set_width_chars(8)
        grid.attach(self.entry_maxcyc, 3, 2, 1, 1)

        grid.attach(Gtk.Label(label="Factor:"), 0, 3, 1, 1)
        self.entry_factor = Gtk.Entry()
        self.entry_factor.set_text("1.0")
        self.entry_factor.set_width_chars(8)
        grid.attach(self.entry_factor, 1, 3, 1, 1)

        grid.attach(Gtk.Label(label="xTB path:"), 0, 4, 1, 1)
        self.entry_xtb = Gtk.Entry()
        self.entry_xtb.set_text("xtb")
        self.entry_xtb.set_hexpand(True)
        grid.attach(self.entry_xtb, 1, 4, 3, 1)
        
        try:
            _XTBCommand = "PDYNAMO3_XTBCOMMAND"
            command = os.getenv ( _XTBCommand )
            self.entry_xtb.set_text(command)
        except:
            pass

        # ---- temporary files folder ----
        # where each fragment's xTB inputs/outputs are written.
        # Empty => use a system temporary folder (tempfile).
        grid.attach(Gtk.Label(label="Temp folder:"), 0, 5, 1, 1)
        self.entry_tmpdir = Gtk.Entry()
        self.entry_tmpdir.set_placeholder_text("(default: system temp)")
        self.entry_tmpdir.set_hexpand(True)
        grid.attach(self.entry_tmpdir, 1, 5, 2, 1)

        self.btn_tmpdir = Gtk.Button(label="Browse...")
        self.btn_tmpdir.connect("clicked", self._on_choose_tmpdir)
        grid.attach(self.btn_tmpdir, 3, 5, 1, 1)

        vbox.pack_start(grid, False, False, 0)

        # ---- progress bar ----
        self.progress = Gtk.ProgressBar()
        self.progress.set_show_text(True)
        self.progress.set_text("idle")
        vbox.pack_start(self.progress, False, False, 0)

        # ---- buttons ----
        row_btn = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        self.btn_run = Gtk.Button(label="Run")
        self.btn_run.connect("clicked", self.on_run)
        self.btn_run.set_sensitive(False)
        row_btn.pack_start(self.btn_run, True, True, 0)

        self.btn_apply = Gtk.Button(label="Apply to system")
        self.btn_apply.connect("clicked", self.on_apply)
        self.btn_apply.set_sensitive(False)
        row_btn.pack_start(self.btn_apply, True, True, 0)

        self.btn_close = Gtk.Button(label="Close")
        self.btn_close.connect("clicked", lambda *_: self.window.destroy())
        row_btn.pack_start(self.btn_close, True, True, 0)

        vbox.pack_start(row_btn, False, False, 0)

        self.window.show_all()
        self._on_level_changed(self.combo_level)  # ajusta visibilidade inicial

    def _on_level_changed(self, combo):
        """Show 'Build fragments' for residue/segment/chain and
        'Import selection fragment' for 'selection' mode."""
        level = combo.get_active_text()
        is_sel = (level == "selection")
        # in selection mode: hide Build, show Import
        self.btn_build.set_visible(not is_sel)
        self.btn_import_sel.set_visible(is_sel)

    def _add_text_column(self, title, col, editable=False):
        renderer = Gtk.CellRendererText()
        if editable:
            renderer.set_property("editable", True)
            renderer.connect("edited", self._on_cell_edited, col)
        column = Gtk.TreeViewColumn(title, renderer, text=col)
        column.set_resizable(True)
        self.treeview.append_column(column)

    # ------------------------------------------------------------------ #
    #  Callbacks                                                          #
    # ------------------------------------------------------------------ #
    def on_build_fragments(self, _button):
        level = self.combo_level.get_active_text()
        if level == "selection":
            # in selection mode fragments come from "Import selection fragment"
            return
        try:
            self.fragments = xfc.build_fragments(self.system, level=level)
        except Exception as e:
            self._error("Failed to build fragments: {}".format(e))
            return

        self.store.clear()
        for frag in self.fragments:
            swatch = make_color_swatch(clustal_rgb(frag.key))
            self.store.append([swatch, True, frag.key, len(frag.atom_indexes),
                               frag.formal_charge, frag.multiplicity, ""])
        self.btn_run.set_sensitive(len(self.fragments) > 0)
        self.btn_apply.set_sensitive(False)
        self.progress.set_text("{} fragments".format(len(self.fragments)))

    def _on_include_toggled(self, _renderer, path):
        """Toggle include/ignore the fragment in the calculation."""
        self.store[path][COL_INCLUDE] = not self.store[path][COL_INCLUDE]

    # ------------------------------------------------------------------ #
    #  Context menu: (de)select all                                       #
    # ------------------------------------------------------------------ #
    def _on_treeview_button_press(self, _widget, event):
        # right button (3) opens the menu
        if event.button == 3:
            menu = Gtk.Menu()
            for label, cb in (("Select all", self._select_all),
                              ("Deselect all", self._deselect_all),
                              ("Invert selection", self._invert_selection)):
                item = Gtk.MenuItem(label=label)
                item.connect("activate", cb)
                menu.append(item)
            menu.show_all()
            menu.popup_at_pointer(event)
            return True
        return False

    def _select_all(self, _item):
        for row in self.store:
            row[COL_INCLUDE] = True

    def _deselect_all(self, _item):
        for row in self.store:
            row[COL_INCLUDE] = False

    def _invert_selection(self, _item):
        for row in self.store:
            row[COL_INCLUDE] = not row[COL_INCLUDE]

    # ------------------------------------------------------------------ #
    #  Import the current vismol selection as a fragment                  #
    # ------------------------------------------------------------------ #
    def on_import_selection(self, _button):
        """Add the current vismol selection as a new fragment in the treeview.

        Can be called multiple times: each call adds a new selection as an
        independent fragment (editable charge/multiplicity).
        """
        indexes = self._get_current_selection_indexes()
        if not indexes:
            self._error("No atoms selected in the viewer.")
            return

        self._selection_counter += 1
        key = "SEL/{}".format(self._selection_counter)
        frag = xfc.Fragment(key, indexes)

        # suggested formal charge = rounded sum of the selection's MM charges
        acc = xfc.SystemAccessor(self.system)
        qsum = sum(acc.mm_charge(i) for i in indexes)
        frag.formal_charge = int(round(qsum))
        frag.multiplicity = 1
        self.fragments.append(frag)

        swatch = make_color_swatch(DEFAULT_RESIDUE_COLOR)  # custom selection: gray
        self.store.append([swatch, True, key, len(indexes),
                           frag.formal_charge, frag.multiplicity, ""])
        self.btn_run.set_sensitive(True)
        self.progress.set_text("{} fragment(s)".format(len(self.fragments)))

    def _get_current_selection_indexes(self):
        """GLOBAL indices (0-based) of the atoms currently selected in vismol.

        NOTE: the vismol atom has a 1-based .index; the engine uses 0-based. Here
        we convert (index - 1). Adjust if your mapping is different.
        """
        try:
            vm = self.easy_session.vm_session
            sel = vm.selections[vm.current_selection]
            indexes = []
            for atom in sel.selected_atoms:
                # atom.index is 1-based in vismol -> 0-based for the engine/pDynamo
                indexes.append(atom.index - 1)
            return sorted(set(indexes))
        except Exception as e:
            self._error("Could not read the current selection: {}".format(e))
            return []

    def _on_cell_edited(self, _renderer, path, new_text, col):
        # validate integer; the treeview is the source of truth before running
        try:
            value = int(new_text)
        except ValueError:
            self._error("Value must be an integer.")
            return
        self.store[path][col] = value

    def on_run(self, _button):
        if self._running:
            return
        # transfer charge/mult/include from the treeview to the fragments (source of truth)
        for i, frag in enumerate(self.fragments):
            frag.formal_charge = int(self.store[i][COL_CHG])
            frag.multiplicity = int(self.store[i][COL_MULT])
            frag.include = bool(self.store[i][COL_INCLUDE])

        if not any(f.include for f in self.fragments):
            self._error("No fragments selected for calculation.")
            return

        try:
            tmpdir = self.entry_tmpdir.get_text().strip()
            params = {
                "method": self.combo_method.get_active_text(),
                "charge_model": self.combo_charge.get_active_text(),
                "boundary_mode": self.combo_boundary.get_active_text(),
                "tolerance": float(self.entry_tol.get_text()),
                "nprocs": int(self.entry_nproc.get_text()),
                "max_cycles": int(self.entry_maxcyc.get_text()),
                "xtb_path": self.entry_xtb.get_text().strip(),
                # temp folder: None => the engine uses the default tempfile
                "workroot": tmpdir if tmpdir else None,
            }
        except ValueError as e:
            self._error("Invalid option: {}".format(e))
            return

        # validate the temp folder, if provided
        if params["workroot"]:
            if not os.path.isdir(params["workroot"]):
                try:
                    os.makedirs(params["workroot"], exist_ok=True)
                except Exception as e:
                    self._error("Cannot create temp folder: {}".format(e))
                    return

        self._running = True
        self.btn_run.set_sensitive(False)
        self.btn_build.set_sensitive(False)
        self.btn_apply.set_sensitive(False)
        self.progress.set_fraction(0.0)
        self.progress.set_text("running...")

        thread = threading.Thread(target=self._run_worker, args=(params,), daemon=True)
        thread.start()

    def _run_worker(self, params):
        """Run the engine in a separate thread; update the UI via GLib.idle_add."""
        def cycle_cb(cycle, max_dq):
            frac = min(cycle / float(params["max_cycles"]), 1.0)
            GLib.idle_add(self._update_progress, frac,
                          "cycle {}: max |dq| = {:.5f}".format(cycle, max_dq))

        try:
            result = xfc.run_self_consistent_parallel(
                self.system, self.fragments, params["xtb_path"],
                method=params["method"], charge_model=params["charge_model"],
                boundary_mode=params["boundary_mode"],
                tolerance=params["tolerance"], max_cycles=params["max_cycles"],
                nprocs=params["nprocs"], workroot=params["workroot"],
                verbose=False, cycle_cb=cycle_cb)
            GLib.idle_add(self._run_done, result, None)
        except Exception as e:
            GLib.idle_add(self._run_done, None, str(e))

    def _update_progress(self, fraction, text):
        self.progress.set_fraction(fraction)
        self.progress.set_text(text)
        return False

    def _run_done(self, result, error):
        self._running = False
        self.btn_build.set_sensitive(True)
        self.btn_run.set_sensitive(True)
        if error:
            self.progress.set_text("error")
            self._error("Run failed: {}".format(error))
            return
        self.result = result
        conv = "converged" if result["converged"] else "NOT converged"
        self.progress.set_fraction(1.0)
        self.progress.set_text("{} in {} cycles (max|dq|={:.5f})".format(
            conv, result["cycles"], result["max_dq"]))
        self.btn_apply.set_sensitive(True)
        return False

    def on_apply(self, _button):
        
        try:
            factor = float(self.entry_factor.get_text())
        except:
            print('Invalid factor value! Using factor = 1.0.')
            factor = 1.0
        
        
        if not self.result:
            return
        try:
            xfc.apply_charges_to_system(self.system, self.result["charges_by_index"], factor)
            total = sum(self.result["charges_by_index"].values())
            self.progress.set_text("applied (total charge = {:.4f})".format(total))
            
        except Exception as e:
            self._error("Failed to apply charges: {}".format(e))

    def _on_destroy(self, *_):
        self.window = None

    def _error(self, message):
        dialog = Gtk.MessageDialog(
            transient_for=self.window, modal=True,
            message_type=Gtk.MessageType.ERROR, buttons=Gtk.ButtonsType.OK,
            text=message)
        dialog.run()
        dialog.destroy()
