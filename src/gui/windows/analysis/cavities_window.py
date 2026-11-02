"""
cavities_window
===============

Analysis > Cavities and Pockets... (2026-10-30): finds the cavities of the
structure shown (util/cavities.py, the two-probe grid method of KVFinder),
lists them with volume, area, depth and lining residues, shows them in the
3D view (surface or grid points, one switch per cavity, in a "Cavities"
object under the system), exports them, and follows the volume of one
cavity along the frames of an object (MD, scan, NEB, umbrella sampling).
"""

import os
import threading
import time
import traceback
import types

import numpy as np

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

from util import cavities as CV
from util.analysis_plot import AnalysisPlot
from gui.widgets.custom_widgets import SystemComboBox, CoordinatesComboBox


def _rgb_hex(rgb):
    return "#{:02x}{:02x}{:02x}".format(*(int(round(255 * c)) for c in rgb))


def add_cavity_object(main, parent, frame, meshes, name):
    """ A vobject under `parent` holding the meshes {cavity id: (surface, points)}
        (util.cavities.cavity_mesh / points_mesh) as representations "cav<id>" and
        "pts<id>", all off; show them with show_cavities(). One frame: `frame` of
        parent. Also used by the DYFFDock window (pockets of the receptor). """
    from vismol.core.vismol_object import VismolObject
    from vismol.libgl.representations import SurfaceRepresentation
    vm_session = main.vm_session
    obj = VismolObject(name=name, index=-1, vismol_session=vm_session, trajectory=[], bonds_pair_of_indexes=[0, 1])
    obj.model_mat = parent.model_mat
    obj.trans_mat = parent.trans_mat
    surfaces = {}
    for cid, (surface, points) in meshes.items():
        if surface is not None:
            surfaces["cav{}".format(cid)] = surface
        if points is not None:
            surfaces["pts{}".format(cid)] = points
    obj.surface_trajectory = [surfaces]
    for key in surfaces:
        rep = SurfaceRepresentation(vismol_object=obj, vismol_glcore=vm_session.vm_glcore, name="surface",
                                    active=False, indexes=[], is_dynamic=False, surface_name=key)
        rep.set_shading_mode("smooth")
        obj.representations[key] = rep
    # . one frame: the cavities of the analysed frame, whatever frame is shown
    obj.frames = np.asarray(parent.frames)[frame:frame + 1]
    obj.active = True
    obj.is_surface = True
    obj.surface_type = "cavities"
    obj.e_id = parent.e_id
    vm_session._add_vismol_object(obj, show_molecule=False, autocenter=False)
    main.main_treeview.add_vismol_object_to_treeview(obj, parent.e_treeview_iter)
    main.add_vobject_to_vobject_liststore_dict(obj)
    return obj


def show_cavities(main, obj, shown, points=False, alpha=0.6):
    """ Turns on the cavities `shown` (ids) of an add_cavity_object() object, as
        surfaces or as grid points. """
    for key, rep in obj.representations.items():
        if not (key.startswith("cav") or key.startswith("pts")):
            continue
        rep.active = int(key[3:]) in shown and (key.startswith("pts") == points)
        rep.set_alpha(1.0 if key.startswith("pts") else alpha)
    main.vm_session.vm_glcore.queue_draw()


def hide_cavity_object(main, obj):
    """ The object stays in the treeview, hidden (a later detection replaces it). """
    if obj is None or obj not in main.vm_session.vm_objects_dic.values():
        return
    obj.active = False
    try:
        main.main_treeview.treestore.set_value(obj.e_treeview_iter, 6, False)
    except Exception:
        pass
    main.vm_session.vm_glcore.queue_draw()


class CavitiesWindow:

    COLUMNS = ("show", "id", "volume", "area", "depth", "buried", "n_res", "residues", "color")

    def __init__(self, main):
        self.main = main
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.window = None
        self.Visible = False
        self._thread = None
        self.result = None          # . last detection (util.cavities result + context)
        self.cavity_object = None   # . the "Cavities" object of the last detection
        self.trajectory = None      # . last volume-per-frame run

    # ------------------------------------------------------------------ window
    def open_window(self):
        if self.window is not None:
            self.window.present()
            return
        self.window = Gtk.Window(title="Cavities and Pockets")
        self.window.set_default_size(900, 680)
        self.window.set_transient_for(self.main.window)
        self.window.connect("destroy", self._on_destroy)
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8, margin=10)
        self.window.add(box)
        box.pack_start(self._build_input(), False, False, 0)
        box.pack_start(self._build_parameters(), False, False, 0)

        run = Gtk.Box(spacing=8)
        self.btn_detect = Gtk.Button(label="Detect Cavities")
        self.btn_detect.connect("clicked", lambda *_: self._safe(self.detect))
        self.label_status = Gtk.Label(xalign=0)
        self.label_status.set_ellipsize(3)
        run.pack_start(self.btn_detect, False, False, 0)
        run.pack_start(self.label_status, True, True, 0)
        box.pack_start(run, False, False, 0)

        self.notebook = Gtk.Notebook()
        self.notebook.append_page(self._build_table_page(), Gtk.Label(label="Cavities"))
        self.notebook.append_page(self._build_trajectory_page(), Gtk.Label(label="Volume along frames"))
        box.pack_start(self.notebook, True, True, 0)

        bottom = Gtk.Box(spacing=8)
        close = Gtk.Button(label="Close")
        close.connect("clicked", lambda *_: self.window.destroy())
        bottom.pack_end(close, False, False, 0)
        box.pack_start(bottom, False, False, 0)
        self.window.show_all()
        self.Visible = True
        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)

    def _build_input(self):
        grid = Gtk.Grid(column_spacing=8, row_spacing=6)
        self.coordinates_combobox = CoordinatesComboBox()
        self.coordinates_combobox.connect("changed", self._on_object_changed)
        self.combobox_systems = SystemComboBox(self.main, coord_combobox=self.coordinates_combobox)
        self.spin_frame = Gtk.SpinButton.new_with_range(0, 0, 1)
        self.spin_frame.set_tooltip_text("The frame analysed by Detect Cavities.")
        grid.attach(Gtk.Label(label="System", xalign=0), 0, 0, 1, 1)
        grid.attach(self.combobox_systems, 1, 0, 1, 1)
        grid.attach(Gtk.Label(label="Object", xalign=0), 2, 0, 1, 1)
        grid.attach(self.coordinates_combobox, 3, 0, 1, 1)
        grid.attach(Gtk.Label(label="Frame", xalign=0), 4, 0, 1, 1)
        grid.attach(self.spin_frame, 5, 0, 1, 1)
        return grid

    def _spin(self, lo, hi, step, value, digits, tip):
        spin = Gtk.SpinButton.new_with_range(lo, hi, step)
        spin.set_digits(digits)
        spin.set_value(value)
        spin.set_tooltip_text(tip)
        return spin

    def _build_parameters(self):
        frame = Gtk.Frame(label=" Parameters ")
        grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
        frame.add(grid)
        self.spin_spacing = self._spin(0.3, 1.5, 0.1, CV.DEFAULT_SPACING, 2,
                                       "Grid step (A). Smaller is more precise and slower (time ~ 1/step^3).")
        self.spin_probe_in = self._spin(0.5, 3.0, 0.1, CV.DEFAULT_PROBE_IN, 2,
                                        "Probe In (A): the space a probe of this radius reaches is the cavity space (1.4 = water).")
        self.spin_probe_out = self._spin(2.0, 15.0, 0.5, CV.DEFAULT_PROBE_OUT, 1,
                                         "Probe Out (A): a large probe rolled from outside defines the bulk solvent. "
                                         "Larger values close wider pockets and find shallower ones.")
        self.spin_removal = self._spin(0.0, 6.0, 0.2, CV.DEFAULT_REMOVAL_DISTANCE, 1,
                                       "Removal distance (A): the layer next to the bulk that does not count as cavity.")
        self.spin_min_volume = self._spin(0.0, 500.0, 1.0, CV.DEFAULT_MIN_VOLUME, 1,
                                          "Cavities smaller than this volume (A^3) are not listed.")
        for k, (text, w) in enumerate((("Grid step (Å)", self.spin_spacing), ("Probe In (Å)", self.spin_probe_in),
                                       ("Probe Out (Å)", self.spin_probe_out), ("Removal distance (Å)", self.spin_removal),
                                       ("Minimum volume (Å³)", self.spin_min_volume))):
            grid.attach(Gtk.Label(label=text, xalign=0), 2 * (k % 3), k // 3, 1, 1)
            grid.attach(w, 2 * (k % 3) + 1, k // 3, 1, 1)

        ignore = Gtk.Box(spacing=10)
        ignore.pack_start(Gtk.Label(label="Leave out:"), False, False, 0)
        self.check_water = Gtk.CheckButton(label="water", active=True)
        self.check_ions = Gtk.CheckButton(label="ions", active=False)
        self.check_ions.set_tooltip_text("Residues of one atom that are not part of the protein (Na+, Cl-, and also metal ions of the active site).")
        self.check_ligands = Gtk.CheckButton(label="ligands", active=False)
        self.check_ligands.set_tooltip_text("Other residues that are neither protein nor water: leave them out to see the empty binding pocket.")
        self.check_hydrogens = Gtk.CheckButton(label="hydrogens", active=False)
        self.check_hydrogens.set_tooltip_text("Crystal structures have none; leave them out to compare a protonated system with one.")
        for w in (self.check_water, self.check_ions, self.check_ligands, self.check_hydrogens):
            ignore.pack_start(w, False, False, 0)
        grid.attach(ignore, 0, 2, 6, 1)

        region = Gtk.Box(spacing=8)
        region.pack_start(Gtk.Label(label="Region:"), False, False, 0)
        self.radio_whole = Gtk.RadioButton.new_with_label(None, "whole molecule")
        self.radio_around = Gtk.RadioButton.new_with_label_from_widget(self.radio_whole, "box around the selected atoms, margin (Å)")
        self.radio_around.set_tooltip_text("Only cavity points inside this box are kept (e.g. select a ligand or the active-site residues). "
                                           "Along frames the box follows those atoms.")
        self.spin_margin = self._spin(0.0, 20.0, 0.5, 4.0, 1, "Margin added to the box of the selected atoms (A).")
        region.pack_start(self.radio_whole, False, False, 0)
        region.pack_start(self.radio_around, False, False, 0)
        region.pack_start(self.spin_margin, False, False, 0)
        grid.attach(region, 0, 3, 6, 1)
        return frame

    def _build_table_page(self):
        page = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=6)
        self.store = Gtk.ListStore(bool, int, float, str, str, str, int, str, str)
        self.view = Gtk.TreeView(model=self.store)
        toggle = Gtk.CellRendererToggle()
        toggle.connect("toggled", self._on_show_toggled)
        self.view.append_column(Gtk.TreeViewColumn("Show", toggle, active=0))
        swatch = Gtk.CellRendererText()
        col = Gtk.TreeViewColumn("#", swatch, text=1, background=8)
        self.view.append_column(col)
        for title, index, fmt in (("Volume (Å³)", 2, "{:.1f}"),):
            cell = Gtk.CellRendererText(xalign=1.0)
            column = Gtk.TreeViewColumn(title, cell)
            column.set_cell_data_func(cell, lambda c, cell, m, it, data: cell.set_property("text", data[1].format(m[it][data[0]])), (index, fmt))
            column.set_sort_column_id(index)
            self.view.append_column(column)
        for title, index in (("Area (Å²)", 3), ("Depth max (Å)", 4), ("Buried", 5), ("Residues", 6)):
            cell = Gtk.CellRendererText(xalign=1.0 if index != 5 else 0.5)
            column = Gtk.TreeViewColumn(title, cell, text=index)
            self.view.append_column(column)
        cell = Gtk.CellRendererText()
        cell.set_property("ellipsize", 3)
        column = Gtk.TreeViewColumn("Lining residues", cell, text=7)
        column.set_expand(True)
        self.view.append_column(column)
        self.view.connect("row-activated", lambda *_: self._safe(self.center_selected))
        self.view.set_tooltip_text("Double-click a cavity to center the view on it.")
        sw = Gtk.ScrolledWindow()
        sw.add(self.view)
        page.pack_start(sw, True, True, 0)

        actions = Gtk.Box(spacing=6)
        b = Gtk.Button(label="Select Lining Residues")
        b.connect("clicked", lambda *_: self._safe(self.select_lining))
        actions.pack_start(b, False, False, 0)
        b = Gtk.Button(label="Center")
        b.connect("clicked", lambda *_: self._safe(self.center_selected))
        actions.pack_start(b, False, False, 0)
        actions.pack_start(Gtk.Label(label="  Show as"), False, False, 0)
        self.combo_display = Gtk.ComboBoxText()
        for text in ("surface", "grid points"):
            self.combo_display.append_text(text)
        self.combo_display.set_active(0)
        self.combo_display.connect("changed", lambda *_: self._apply_display())
        actions.pack_start(self.combo_display, False, False, 0)
        actions.pack_start(Gtk.Label(label="opacity"), False, False, 0)
        self.spin_alpha = self._spin(0.1, 1.0, 0.1, 0.6, 1, "Opacity of the cavity surfaces.")
        self.spin_alpha.connect("value-changed", lambda *_: self._apply_display())
        actions.pack_start(self.spin_alpha, False, False, 0)
        b = Gtk.Button(label="Export...")
        b.set_tooltip_text("The table (CSV) and the cavity points (PDB, residue CAV, residue number = cavity).")
        b.connect("clicked", lambda *_: self._safe(self.export))
        actions.pack_end(b, False, False, 0)
        page.pack_start(actions, False, False, 0)
        return page

    def _build_trajectory_page(self):
        page = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=6)
        grid = Gtk.Grid(column_spacing=8, row_spacing=6)
        self.spin_first = self._spin(0, 0, 1, 0, 0, "First frame.")
        self.spin_last = self._spin(0, 0, 1, 0, 0, "Last frame.")
        self.spin_stride = self._spin(1, 100000, 1, 1, 0, "Use every n-th frame.")
        import multiprocessing
        self.spin_procs = self._spin(1, multiprocessing.cpu_count(), 1, max(1, multiprocessing.cpu_count() - 1), 0,
                                     "Frames analysed at the same time (processes).")
        self.spin_overlap = self._spin(0.05, 1.0, 0.05, 0.3, 2,
                                       "A cavity of a frame is the followed pocket when this fraction of its lining "
                                       "atoms lines the pocket in the detection frame.")
        for k, (text, w) in enumerate((("First", self.spin_first), ("Last", self.spin_last), ("Stride", self.spin_stride),
                                       ("Processes", self.spin_procs), ("Min. overlap", self.spin_overlap))):
            grid.attach(Gtk.Label(label=text, xalign=0), 2 * k, 0, 1, 1)
            grid.attach(w, 2 * k + 1, 0, 1, 1)
        self.btn_follow = Gtk.Button(label="Follow Selected Cavity")
        self.btn_follow.set_tooltip_text("Volume of the cavity selected in the table in each frame of the object, "
                                         "with the parameters above. The pocket is recognised by the atoms lining it.")
        self.btn_follow.connect("clicked", lambda *_: self._safe(self.follow))
        grid.attach(self.btn_follow, 10, 0, 1, 1)
        page.pack_start(grid, False, False, 0)
        self.plot = AnalysisPlot()
        self.plot.set_size_request(-1, 260)
        page.pack_start(self.plot, True, True, 0)
        hb = Gtk.Box(spacing=6)
        self.label_traj = Gtk.Label(xalign=0)
        self.label_traj.set_line_wrap(True)
        hb.pack_start(self.label_traj, True, True, 0)
        b = Gtk.Button(label="Save Data...")
        b.connect("clicked", lambda *_: self._safe(self.save_trajectory))
        hb.pack_end(b, False, False, 0)
        page.pack_start(hb, False, False, 0)
        return page

    # ------------------------------------------------------------------ input
    def _on_object_changed(self, *_):
        vobject = self._vobject(quiet=True)
        n = len(vobject.frames) if vobject is not None else 1
        for spin in (self.spin_frame, self.spin_first, self.spin_last):
            spin.set_range(0, max(0, n - 1))
        try:
            self.spin_frame.set_value(min(int(self.vm_session.get_frame() or 0), n - 1))
        except Exception:
            self.spin_frame.set_value(0)
        self.spin_last.set_value(max(0, n - 1))

    def _vobject(self, quiet=False):
        try:
            return self.vm_session.vm_objects_dic[self.coordinates_combobox.get_vobject_id()]
        except Exception:
            if quiet:
                return None
            raise ValueError("Choose a system and an object.")

    def _atoms(self, vobject):
        """ Indices of the atoms used (the options "leave out"), their symbols and residue labels. """
        keep, symbols, labels = [], [], []
        for i in sorted(vobject.atoms):
            atom = vobject.atoms[i]
            residue = atom.residue
            symbol = (atom.symbol or "").capitalize()
            labels.append((getattr(atom.chain, "name", "") if atom.chain is not None else "",
                           getattr(residue, "name", "") if residue is not None else "",
                           getattr(residue, "index", 0) if residue is not None else 0))
            symbols.append(symbol)
            if residue is not None:
                if self.check_water.get_active() and residue.is_solvent:
                    continue
                if not residue.is_protein and not residue.is_solvent:
                    single = len(residue.atoms) == 1
                    if single and self.check_ions.get_active():
                        continue
                    if not single and self.check_ligands.get_active():
                        continue
            if self.check_hydrogens.get_active() and symbol == "H":
                continue
            keep.append(i)
        return np.array(keep, dtype=np.int64), symbols, labels

    def _options(self):
        return {"spacing": self.spin_spacing.get_value(), "probe_in": self.spin_probe_in.get_value(),
                "probe_out": self.spin_probe_out.get_value(), "removal_distance": self.spin_removal.get_value(),
                "min_volume": self.spin_min_volume.get_value()}

    def _selected_indices(self, vobject):
        sel = self.vm_session.selections[self.vm_session.current_selection]
        indices = sorted(a.index - 1 for a in sel.selected_atoms if a.vm_object is vobject)
        if not indices:
            raise ValueError("No atoms of this object are selected (the region is 'box around the selected atoms').")
        return indices

    # ------------------------------------------------------------------ detection
    def detect(self):
        if self._thread is not None:
            return
        vobject = self._vobject()
        frame = int(self.spin_frame.get_value())
        keep, symbols, labels = self._atoms(vobject)
        if len(keep) == 0:
            raise ValueError("No atoms left.")
        frames = np.asarray(vobject.frames)
        coords = np.asarray(frames[min(frame, len(frames) - 1)], dtype=np.float64)[keep]
        table = self.vm_session.periodic_table
        radii = CV.radii_for([symbols[i] for i in keep], fallback=lambda s: table.get_vdw_radius(symbol=s))
        options = self._options()
        box_atoms = None
        if self.radio_around.get_active():
            selected = self._selected_indices(vobject)
            box_atoms = (selected, self.spin_margin.get_value())
            options["box"] = CV.box_around(np.asarray(frames[min(frame, len(frames) - 1)])[selected], box_atoms[1])
        self.btn_detect.set_sensitive(False)
        self._t0 = time.time()
        self.label_status.set_text("running...")

        def work():
            try:
                result = CV.detect_cavities(coords, radii, **options)
                meshes = {}
                for c in result["cavities"]:
                    meshes[c["id"]] = (CV.cavity_mesh(result, c["id"]), CV.points_mesh(result, c["id"]))
                GLib.idle_add(self._detected, result, meshes, None)
            except Exception as error:
                traceback.print_exc()
                GLib.idle_add(self._detected, None, None, str(error))

        self._context = {"vobject": vobject, "frame": frame, "keep": keep, "labels": labels, "radii": radii,
                         "options": options, "box_atoms": box_atoms}
        self._thread = threading.Thread(target=work, daemon=True)
        self._thread.start()

    def _detected(self, result, meshes, error):
        self._thread = None
        self.btn_detect.set_sensitive(True)
        if error:
            self.label_status.set_text("failed")
            self._message(error)
            return False
        ctx = self._context
        keep, labels = ctx["keep"], ctx["labels"]
        for c in result["cavities"]:
            c["lining_full"] = keep[c["lining_atoms"]] if len(c["lining_atoms"]) else np.zeros(0, dtype=np.int64)
            c["residues"] = CV.lining_residues({"lining_atoms": c["lining_full"]}, labels)
        result["context"] = ctx
        self.result = result
        self.store.clear()
        for c in result["cavities"]:
            color = CV.CAVITY_COLORS[(c["id"] - 1) % len(CV.CAVITY_COLORS)]
            self.store.append([c["id"] <= 10, c["id"], c["volume"],
                               "" if c["area"] is None else "{:.1f}".format(c["area"]),
                               "-" if c["buried"] else "{:.1f}".format(c["depth_max"]),
                               "yes" if c["buried"] else "no", len(c["residues"]),
                               " ".join("{}{}".format(name, num) for _ch, name, num in c["residues"]), _rgb_hex(color)])
        total = sum(c["volume"] for c in result["cavities"])
        self.label_status.set_text("{} cavities, {:.0f} Å³ in all, frame {} of {} ({} atoms), {:.1f} s".format(
            len(result["cavities"]), total, ctx["frame"], ctx["vobject"].name, len(keep), time.time() - self._t0))
        if result["cavities"]:
            self._make_object(meshes)
            self.view.get_selection().select_path(Gtk.TreePath(0))
        return False

    # ------------------------------------------------------------------ 3D view
    def _make_object(self, meshes):
        # . the previous detection stays in the treeview, hidden
        hide_cavity_object(self.main, self.cavity_object)
        self.cavity_object = add_cavity_object(self.main, self._context["vobject"], self._context["frame"], meshes,
                                               "Cavities (frame {})".format(self._context["frame"]))
        self._apply_display()

    def _apply_display(self):
        obj = self.cavity_object
        if obj is None or self.result is None:
            return
        show_cavities(self.main, obj, {row[1] for row in self.store if row[0]},
                      points=self.combo_display.get_active() == 1, alpha=self.spin_alpha.get_value())

    def _on_show_toggled(self, _cell, path):
        self.store[path][0] = not self.store[path][0]
        self._apply_display()

    def _selected_cavity(self):
        if self.result is None:
            raise ValueError("Detect the cavities first.")
        model, it = self.view.get_selection().get_selected()
        if it is None:
            raise ValueError("Select a cavity in the table.")
        return self.result["cavities"][model[it][1] - 1]

    def select_lining(self):
        c = self._selected_cavity()
        vobject = self.result["context"]["vobject"]
        atoms = set()
        residues = {(ch, name, num) for ch, name, num in c["residues"]}
        for i, atom in vobject.atoms.items():
            residue = atom.residue
            key = (getattr(atom.chain, "name", "") if atom.chain is not None else "",
                   getattr(residue, "name", "") if residue is not None else "",
                   getattr(residue, "index", 0) if residue is not None else 0)
            if key in residues:
                atoms.add(i)
        selection = self.vm_session.selections[self.vm_session.current_selection]
        # . selecting_by_indexes(clear=True) keeps the old atoms while the selection is
        #   active (vismol _clear_selection_buffer): clear the flags here
        for vm_object in self.vm_session.vm_objects_dic.values():
            for atom in vm_object.atoms.values():
                atom.selected = False
        selection.selecting_by_indexes(vismol_object=vobject, indexes=sorted(atoms), clear=False)
        selection.active = True
        self.vm_session.vm_glcore.queue_draw()

    def center_selected(self):
        c = self._selected_cavity()
        vobject = self.result["context"]["vobject"]
        self.vm_session.vm_glcore.center_on_coordinates(vobject, np.asarray(c["centre"], dtype=np.float32))

    # ------------------------------------------------------------------ files
    def export(self):
        if self.result is None:
            raise ValueError("Detect the cavities first.")
        dialog = Gtk.FileChooserDialog(title="Export cavities (base name)", parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_current_name("cavities_{}".format(self.result["context"]["vobject"].name.replace(" ", "_")))
        path = dialog.get_filename() if dialog.run() == Gtk.ResponseType.OK else None
        dialog.destroy()
        if not path:
            return
        base = os.path.splitext(path)[0]
        self.export_files(base)
        self._message("Written:\n{0}.csv (table)\n{0}.pdb (cavity points)".format(base))

    def export_files(self, base):
        labels = self.result["context"]["labels"]
        lining = {c["id"]: c["lining_atoms"] for c in self.result["cavities"]}
        for c in self.result["cavities"]:
            c["lining_atoms"] = c["lining_full"]
        try:
            CV.write_table_csv(self.result, base + ".csv", atom_labels=labels)
        finally:
            for c in self.result["cavities"]:
                c["lining_atoms"] = lining[c["id"]]
        CV.write_points_pdb(self.result, base + ".pdb")
        return base + ".csv", base + ".pdb"

    # ------------------------------------------------------------------ trajectory
    def follow(self):
        if self._thread is not None:
            return
        c = self._selected_cavity()
        ctx = self.result["context"]
        vobject = ctx["vobject"]
        first, last, stride = int(self.spin_first.get_value()), int(self.spin_last.get_value()), int(self.spin_stride.get_value())
        frame_ids = list(range(first, last + 1, max(1, stride)))
        if len(frame_ids) < 2:
            raise ValueError("Choose at least two frames (First, Last, Stride).")
        all_frames = np.asarray(vobject.frames)
        frames = [np.asarray(all_frames[f], dtype=np.float64)[ctx["keep"]] for f in frame_ids]
        options = dict(ctx["options"])
        options.pop("box", None)
        if ctx["box_atoms"] is not None:
            # . the selected atoms, as positions among the atoms used
            position = {int(a): k for k, a in enumerate(ctx["keep"])}
            atoms = [position[a] for a in ctx["box_atoms"][0] if a in position]
            options["box_atoms"] = (atoms, ctx["box_atoms"][1])
        reference = c["lining_atoms"]
        processes = int(self.spin_procs.get_value())
        overlap = self.spin_overlap.get_value()
        self.btn_follow.set_sensitive(False)
        self._t0 = time.time()
        self.label_traj.set_text("{} frames, running...".format(len(frame_ids)))

        def work():
            try:
                rows = CV.pocket_volume_trajectory(frames, ctx["radii"], reference, options=options,
                                                   min_overlap=overlap, processes=processes)
                GLib.idle_add(self._followed, c["id"], frame_ids, rows, None)
            except Exception as error:
                traceback.print_exc()
                GLib.idle_add(self._followed, c["id"], frame_ids, None, str(error))

        self._thread = threading.Thread(target=work, daemon=True)
        self._thread.start()

    def _followed(self, cavity_id, frame_ids, rows, error):
        self._thread = None
        self.btn_follow.set_sensitive(True)
        if error:
            self.label_traj.set_text("failed")
            self._message(error)
            return False
        volumes = [r["volume"] for r in rows]
        self.trajectory = {"cavity": cavity_id, "frames": frame_ids, "rows": rows,
                           "object": self.result["context"]["vobject"].name}
        self.plot.set_lines([{"x": frame_ids, "y": volumes, "label": "cavity {}".format(cavity_id),
                              "color": CV.CAVITY_COLORS[(cavity_id - 1) % len(CV.CAVITY_COLORS)]}],
                            xlabel="Frame", ylabel="Volume (Å³)", title="Cavity {} along {}".format(
                                cavity_id, self.result["context"]["vobject"].name), ymin=0.0)
        lost = sum(1 for r in rows if not r["cavities"])
        v = np.array(volumes)
        text = "Cavity {}: {:.0f} ± {:.0f} Å³ (min {:.0f}, max {:.0f}) over {} frames, {:.1f} s.".format(
            cavity_id, v.mean(), v.std(), v.min(), v.max(), len(v), time.time() - self._t0)
        if lost:
            text += " In {} frames no cavity lines the pocket (volume 0: closed, or merged below the overlap).".format(lost)
        self.label_traj.set_text(text)
        return False

    def save_trajectory(self):
        if self.trajectory is None:
            raise ValueError("Follow a cavity first.")
        dialog = Gtk.FileChooserDialog(title="Save data", parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name("cavity{}_volume.dat".format(self.trajectory["cavity"]))
        path = dialog.get_filename() if dialog.run() == Gtk.ResponseType.OK else None
        dialog.destroy()
        if path:
            self.write_trajectory(path)

    def write_trajectory(self, path):
        t = self.trajectory
        with open(path, "w") as f:
            f.write("# EasyHybrid Cavities: volume of cavity {} along {}\n".format(t["cavity"], t["object"]))
            f.write("# frame  volume_A3  matched_cavities  overlap  all_cavities_volume_A3\n")
            for frame, r in zip(t["frames"], t["rows"]):
                f.write("{:6d} {:10.1f}  {:<16s} {:7.2f} {:10.1f}\n".format(
                    frame, r["volume"], ",".join(str(i) for i in r["cavities"]) or "-", r["overlap"], r["total_volume"]))
        return path

    # ------------------------------------------------------------------ misc
    def update(self, *args, **kwargs):
        """ Called by the main window when systems or objects change. """
        if not self.Visible:
            return
        sid = self.combobox_systems.get_system_id()
        if sid is not None and sid in self.main.vobject_liststore_dict:
            self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[sid])

    def _safe(self, callback):
        try:
            callback()
        except Exception as error:
            traceback.print_exc()
            self._message(str(error))

    def _message(self, text):
        from automation import events as automation_events
        if automation_events.active():
            automation_events.record("info", text, title="Cavities")
            return
        dialog = Gtk.MessageDialog(transient_for=self.window, modal=True, message_type=Gtk.MessageType.INFO,
                                   buttons=Gtk.ButtonsType.OK, text=text)
        dialog.run()
        dialog.destroy()

    def _on_destroy(self, *_):
        self.window = None
        self.Visible = False
