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

Janela GTK (Fase 2) para a ferramenta de cargas MM auto-consistentes por
fragmento (motor em util/xtb_fragment_charges.py).

Fluxo de uso:
  1. escolher o nivel (residuo/segmento/cadeia) -> "Build fragments" popula a
     treeview com carga formal e multiplicidade SUGERIDAS;
  2. o usuario pode editar carga/multiplicidade de cada fragmento na treeview
     (a treeview e' a fonte da verdade antes de rodar);
  3. escolher metodo (GFN0/1/2), tipo de carga (CM5/Mulliken), tolerancia,
     modo de fronteira, path do xTB e nº de processos;
  4. "Run" -> roda o loop auto-consistente em THREAD separada (cada fragmento
     e' um processo, via multiprocessing.Pool); a barra de progresso atualiza
     por ciclo;
  5. "Apply to system" -> grava as cargas finais em system.mmState.charges.

Integracao: instanciar com (main, system) e chamar open_window(). Ex.:
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
#  Paleta de cores Clustal X (por tipo de aminoacido), em RGB 0-255.           #
#  Versao simplificada "por residuo" (sem as regras dependentes de coluna de   #
#  alinhamento): cada grupo de residuos recebe uma cor fixa.                    #
#  [VERIFY: confira se estas sao as cores Clustal que voce quer no paper]       #
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
DEFAULT_RESIDUE_COLOR = (200, 200, 200)  # cinza para o que nao casar

# codigo 3-letras -> 1-letra (para casar com a paleta a partir do nome do residuo)
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
    # variantes de protonacao comuns -> mesmo residuo base
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "CYM": "C",
    "ASH": "D", "GLH": "E", "LYN": "K",
}


def _residue_one_letter(fragment_key):
    """Extrai o codigo de 1 letra do residuo a partir da chave 'CHAIN/RESNAME/SEQ'."""
    try:
        resname = fragment_key.split("/")[1].upper()
    except Exception:
        return None
    if len(resname) == 1:
        return resname
    return _THREE_TO_ONE.get(resname[:3])


def clustal_rgb(fragment_key):
    """Cor RGB (0-255) Clustal para a chave de fragmento; cinza se desconhecido."""
    one = _residue_one_letter(fragment_key)
    return CLUSTAL_COLORS.get(one, DEFAULT_RESIDUE_COLOR)


def make_color_swatch(rgb, width=18, height=18):
    """Cria um GdkPixbuf.Pixbuf preenchido com a cor rgb (0-255)."""
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
    """Janela da ferramenta de cargas por fragmento (xTB)."""

    def __init__(self, main=None, system=None):
        self.main = main.window
        self.easy_session = main
        self.system = system
        self.window = None
        self.fragments = []          # lista de xfc.Fragment
        self.result = None           # dict devolvido pelo motor
        self._running = False
        self._selection_counter = 0  # numera os fragmentos importados de selecao

    # ------------------------------------------------------------------ #
    #  Construcao da UI                                                    #
    # ------------------------------------------------------------------ #
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

        # ---- linha 1: nivel + build / import selection ----
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

        # botao para importar a selecao atual do vismol como um fragmento
        # (visivel apenas no modo 'selection')
        self.btn_import_sel = Gtk.Button(label="Import selection fragment")
        self.btn_import_sel.connect("clicked", self.on_import_selection)
        row1.pack_start(self.btn_import_sel, False, False, 0)
        vbox.pack_start(row1, False, False, 0)

        # ---- treeview de fragmentos (editavel) ----
        # store: color(pixbuf), include(bool), key, natoms, charge, mult, status
        from gi.repository import GdkPixbuf
        self.store = Gtk.ListStore(GdkPixbuf.Pixbuf, bool, str, int, int, int, str)
        self.treeview = Gtk.TreeView(model=self.store)

        # coluna 0: swatch de cor (Clustal)
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

        # menu de contexto (clique direito) para (de)selecionar todos
        self.treeview.connect("button-press-event", self._on_treeview_button_press)

        # ---- opcoes de calculo ----
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

        vbox.pack_start(grid, False, False, 0)

        # ---- barra de progresso ----
        self.progress = Gtk.ProgressBar()
        self.progress.set_show_text(True)
        self.progress.set_text("idle")
        vbox.pack_start(self.progress, False, False, 0)

        # ---- botoes ----
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
        """Mostra 'Build fragments' para residue/segment/chain e
        'Import selection fragment' para o modo 'selection'."""
        level = combo.get_active_text()
        is_sel = (level == "selection")
        # no modo selection: some o Build, aparece o Import
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
            # no modo selection os fragmentos vem de "Import selection fragment"
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
        """Alterna incluir/ignorar o fragmento no calculo."""
        self.store[path][COL_INCLUDE] = not self.store[path][COL_INCLUDE]

    # ------------------------------------------------------------------ #
    #  Menu de contexto: (de)selecionar todos                             #
    # ------------------------------------------------------------------ #
    def _on_treeview_button_press(self, _widget, event):
        # botao direito (3) abre o menu
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
    #  Importar a selecao atual do vismol como um fragmento               #
    # ------------------------------------------------------------------ #
    def on_import_selection(self, _button):
        """Adiciona a selecao atual do vismol como um novo fragmento na treeview.

        Pode ser chamado varias vezes: cada chamada adiciona uma nova selecao
        como um fragmento independente (carga/multiplicidade editaveis).
        """
        indexes = self._get_current_selection_indexes()
        if not indexes:
            self._error("No atoms selected in the viewer.")
            return

        self._selection_counter += 1
        key = "SEL/{}".format(self._selection_counter)
        frag = xfc.Fragment(key, indexes)

        # carga formal sugerida = soma das cargas MM da selecao, arredondada
        acc = xfc.SystemAccessor(self.system)
        qsum = sum(acc.mm_charge(i) for i in indexes)
        frag.formal_charge = int(round(qsum))
        frag.multiplicity = 1
        self.fragments.append(frag)

        swatch = make_color_swatch(DEFAULT_RESIDUE_COLOR)  # selecao custom: cinza
        self.store.append([swatch, True, key, len(indexes),
                           frag.formal_charge, frag.multiplicity, ""])
        self.btn_run.set_sensitive(True)
        self.progress.set_text("{} fragment(s)".format(len(self.fragments)))

    def _get_current_selection_indexes(self):
        """Indices GLOBAIS (base-0) dos atomos atualmente selecionados no vismol.

        ATENCAO: o atom do vismol tem .index base-1; o motor usa base-0. Aqui
        convertemos (index - 1). Ajuste se o seu mapeamento for diferente.
        """
        try:
            vm = self.easy_session.vm_session
            sel = vm.selections[vm.current_selection]
            indexes = []
            for atom in sel.selected_atoms:
                # atom.index e' base-1 no vismol -> base-0 para o motor/pDynamo
                indexes.append(atom.index - 1)
            return sorted(set(indexes))
        except Exception as e:
            self._error("Could not read the current selection: {}".format(e))
            return []

    def _on_cell_edited(self, _renderer, path, new_text, col):
        # valida inteiro; a treeview e' a fonte da verdade antes de rodar
        try:
            value = int(new_text)
        except ValueError:
            self._error("Value must be an integer.")
            return
        self.store[path][col] = value

    def on_run(self, _button):
        if self._running:
            return
        # transfere carga/mult/include da treeview para os fragmentos (fonte da verdade)
        for i, frag in enumerate(self.fragments):
            frag.formal_charge = int(self.store[i][COL_CHG])
            frag.multiplicity = int(self.store[i][COL_MULT])
            frag.include = bool(self.store[i][COL_INCLUDE])

        if not any(f.include for f in self.fragments):
            self._error("No fragments selected for calculation.")
            return

        try:
            params = {
                "method": self.combo_method.get_active_text(),
                "charge_model": self.combo_charge.get_active_text(),
                "boundary_mode": self.combo_boundary.get_active_text(),
                "tolerance": float(self.entry_tol.get_text()),
                "nprocs": int(self.entry_nproc.get_text()),
                "max_cycles": int(self.entry_maxcyc.get_text()),
                "xtb_path": self.entry_xtb.get_text().strip(),
            }
        except ValueError as e:
            self._error("Invalid option: {}".format(e))
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
        """Roda o motor em thread separada; atualiza a UI via GLib.idle_add."""
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
                nprocs=params["nprocs"], verbose=False, cycle_cb=cycle_cb)
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
