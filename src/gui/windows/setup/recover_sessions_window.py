#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Startup / crash recovery (2026-11-01): a session that is changed but was
#  never given a file name yet autosaves to <easyhybrid_tmp>/<random_code>.easy
#  (see io_data.py's register_change_and_maybe_autosave/_do_autosave/
#  save_easyhybrid_session). Once the user later does Save As, or the app is
#  simply closed/killed without ever saving, that file is left behind -- on a
#  real install, hundreds of them over time (913 MB / 175 files found in this
#  project's own .tmp/ while building this). This window lists whatever is
#  still there at startup, with an on-demand preview (unpickled only when a
#  row is selected, since some of these can be very large) and lets the user
#  open, keep a copy of, delete or ignore them.

import os
import pickle
import shutil
import time

import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk  # noqa: E402

from automation import events as automation_events  # noqa: E402


def _format_size(num_bytes):
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return "{:.0f} {}".format(size, unit) if unit == "B" else "{:.1f} {}".format(size, unit)
        size /= 1024.0


def _fast_preview_unpickler(fileobj):
    """ A pickle.Unpickler that resolves pMolecule.System to a subclass
        whose __setstate__ skips Connectivity.FromAtoms() -- bond-order/
        aromaticity perception from scratch (graph pattern matching over
        every ring), which pMolecule.System.__setstate__ always runs, and
        which turned out to be its ENTIRE unpickling cost, regardless of
        how much trajectory data the file holds: profiled (cProfile) on a
        147 MB .tmp/*.easy file while speeding this preview up -- 1.13 s
        of a 1.21 s total load, against ~0.1 s for the very same file once
        this is skipped. A preview only ever reads label/e_tag/atom
        count/force field/job history, never bonds or rings, so this is
        safe -- but ONLY for a preview: the System objects this produces
        have no working connectivity (._connectivity is None) and must
        never be used for anything else. """
    from pMolecule.System import System as RealSystem
    from pMolecule.Atom import AtomContainer
    from pMolecule.Sequence import Sequence

    class _FastPreviewSystem(RealSystem):
        def __setstate__(self, state):
            if "sequence" in state:
                sequence = Sequence.FromMapping(state.pop("sequence"))
                atoms = AtomContainer.FromIterable(sequence.GatherAtoms(), attributes=state.pop("atoms", None))
            else:
                atoms = AtomContainer.FromMapping(state.pop("atoms", None))
                sequence = None
            state.pop("bonds", None)                        # . never reconstructed -- see docstring above
            self.__dict__["_atoms"] = atoms
            self.__dict__["_connectivity"] = None
            self.__dict__["_sequence"] = sequence
            for key, value in state.items():
                if value is not None:
                    self.__dict__[key] = value

    class _FastPreviewUnpickler(pickle.Unpickler):
        def find_class(self, module, name):
            if module == "pMolecule.System" and name == "System":
                return _FastPreviewSystem
            return super().find_class(module, name)

    return _FastPreviewUnpickler(fileobj)


def describe_tmp_easy_file(path):
    """ Unpickles `path` (an EasyHybrid .easy autosave) and returns a short,
        plain-text description: one block per system (label, tag, atom
        count, force field, its vismol objects) and its job history, if
        any. Never raises -- a corrupted, truncated, or just-too-old-a-
        format file gets a one-line explanation instead of a traceback. """
    try:
        from pdynamo.pDynamo2EasyHybrid.io_data import _ensure_numpy_pickle_compat
        _ensure_numpy_pickle_compat()
        with open(path, "rb") as f:
            data = _fast_preview_unpickler(f).load()
    except Exception as error:
        return "Could not read this file:\n{}".format(error)

    systems = data.get("systems") or []
    if not systems:
        return "No system in this session."

    blocks = []
    for entry in systems:
        system = entry.get("system")
        if system is None:
            continue
        label = getattr(system, "label", "?")
        tag = getattr(system, "e_tag", "")
        lines = ["System: {} ({})".format(label, tag) if tag else "System: {}".format(label)]
        try:
            natoms = len(system.atoms)
        except Exception:
            natoms = "?"
        try:
            force_field = getattr(system.mmModel, "forceField", None)
        except Exception:
            force_field = None
        qc_model = getattr(system, "qcModel", None)
        info = "{} atoms".format(natoms)
        if force_field:
            info += "  -  force field: {}".format(force_field)
        if qc_model is not None:
            info += "  -  QC model: {}".format(type(qc_model).__name__)
        lines.append("  " + info)

        vobjects = entry.get("vobjects") or []
        if vobjects:
            names = ", ".join(v.get("name", "?") for v in vobjects)
            lines.append("  {} object(s): {}".format(len(vobjects), names))

        history = getattr(system, "e_job_history", None) or {}
        if history:
            lines.append("  {} job(s) in the history:".format(len(history)))
            for step in sorted(history):
                job = history[step]
                job_line = "    [{}] {} - {}".format(step, job.get("simulation_type", "?"), job.get("status", "?"))
                if job.get("remote"):
                    job_line += " (remote: {})".format(job["remote"])
                lines.append(job_line)
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks)


class RecoverSessionsWindow:
    def __init__(self, main=None):
        self.main = main
        self.window = None
        self.liststore = None
        self.treeview = None
        self.preview_view = None

    # ------------------------------------------------------------------ scan
    def _tmp_folder(self):
        vm_session = getattr(self.main, "vm_session", None)
        vm_config = getattr(vm_session, "vm_config", None)
        return getattr(vm_config, "easyhybrid_tmp", None)

    def _list_files(self):
        folder = self._tmp_folder()
        if not folder or not os.path.isdir(folder):
            return []
        out = []
        for name in os.listdir(folder):
            if not name.endswith(".easy"):
                continue
            path = os.path.join(folder, name)
            try:
                stat = os.stat(path)
            except OSError:
                continue
            out.append((path, stat.st_mtime, stat.st_size))
        out.sort(key=lambda item: item[1], reverse=True)
        return out

    # ------------------------------------------------------------------ startup entry point
    def maybe_open_at_startup(self):
        """ Called once, right after the main window is up. Opens this
            window only when there is at least one orphaned .tmp/*.easy
            file; in automation mode (control server with dialogs off) it
            never opens a window -- an event is recorded instead, same
            convention as simple_dialog.py's info()/error(). """
        files = self._list_files()
        if not files:
            return
        if automation_events.active():
            automation_events.record(
                "info", "{} unsaved session file(s) in {}".format(len(files), self._tmp_folder()),
                title="Unsaved sessions")
            return
        self.open_window()

    # ------------------------------------------------------------------ build
    def open_window(self):
        if self.window is not None:
            self._reload()
            if self.window is not None:
                self.window.present()
            return
        self.window = Gtk.Window(title="Unsaved Sessions")
        if self.main is not None and getattr(self.main, "window", None) is not None:
            self.window.set_transient_for(self.main.window)
        self.window.set_default_size(880, 480)
        self.window.connect("destroy", self._on_destroy)
        outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8, border_width=10)
        self.window.add(outer)

        label = Gtk.Label(xalign=0)
        label.set_markup(
            "<b>EasyHybrid found session files that were never saved</b> -- autosaves of a "
            "session that had not been given a file name yet. Select one to preview its "
            "contents, then <i>Open</i>, <i>Save As...</i> or <i>Delete</i> it; or clear/ignore "
            "the whole list.")
        label.set_line_wrap(True)
        outer.pack_start(label, False, False, 0)

        paned = Gtk.Paned(orientation=Gtk.Orientation.HORIZONTAL)
        paned.set_wide_handle(True)
        outer.pack_start(paned, True, True, 0)

        self.liststore = Gtk.ListStore(str, str, str, str)      # path, name, modified, size
        self.treeview = Gtk.TreeView(model=self.liststore)
        for index, title in enumerate(("File", "Modified", "Size")):
            renderer = Gtk.CellRendererText()
            column = Gtk.TreeViewColumn(title, renderer, text=index + 1)
            column.set_resizable(True)
            column.set_sort_column_id(index + 1)
            self.treeview.append_column(column)
        self.treeview.get_selection().connect("changed", self._on_selection_changed)
        scroller_left = Gtk.ScrolledWindow(min_content_width=340)
        scroller_left.add(self.treeview)
        paned.pack1(scroller_left, True, False)

        self.preview_view = Gtk.TextView(editable=False, cursor_visible=False, monospace=True,
                                         wrap_mode=Gtk.WrapMode.WORD, left_margin=6, right_margin=6, top_margin=6)
        scroller_right = Gtk.ScrolledWindow()
        scroller_right.add(self.preview_view)
        paned.pack2(scroller_right, True, False)

        buttons = Gtk.Box(spacing=6)
        outer.pack_start(buttons, False, False, 0)
        self.button_open = Gtk.Button(label="Open")
        self.button_open.connect("clicked", self._on_open)
        self.button_save_as = Gtk.Button(label="Save As...")
        self.button_save_as.connect("clicked", self._on_save_as)
        self.button_delete = Gtk.Button(label="Delete")
        self.button_delete.connect("clicked", self._on_delete)
        for button in (self.button_open, self.button_save_as, self.button_delete):
            buttons.pack_start(button, False, False, 0)
        buttons.pack_start(Gtk.Box(hexpand=True), True, True, 0)
        button_clear = Gtk.Button(label="Clear All")
        button_clear.connect("clicked", self._on_clear_all)
        buttons.pack_start(button_clear, False, False, 0)
        button_ignore = Gtk.Button(label="Ignore")
        button_ignore.connect("clicked", lambda *_: self.window.destroy())
        buttons.pack_start(button_ignore, False, False, 0)

        self._reload()
        if self.window is not None:
            self.window.show_all()

    def _reload(self):
        self.liststore.clear()
        files = self._list_files()
        for path, mtime, size in files:
            modified = time.strftime("%d/%m/%Y %H:%M", time.localtime(mtime))
            self.liststore.append([path, os.path.basename(path), modified, _format_size(size)])
        self._set_buttons_sensitive(False)
        self.preview_view.get_buffer().set_text("Select a file to preview its contents.")
        if not files:
            self.window.destroy()               # . nothing left to recover

    def _set_buttons_sensitive(self, flag):
        for button in (self.button_open, self.button_save_as, self.button_delete):
            button.set_sensitive(flag)

    def _selected_path(self):
        model, treeiter = self.treeview.get_selection().get_selected()
        if treeiter is None:
            return None
        return model[treeiter][0]

    # ------------------------------------------------------------------ handlers
    def _on_selection_changed(self, selection):
        path = self._selected_path()
        self._set_buttons_sensitive(path is not None)
        if path is None:
            return
        self.preview_view.get_buffer().set_text(describe_tmp_easy_file(path))

    def _confirm_discard_current_session(self):
        """ The current (pre-Open) session is almost always the fresh, empty
            one from this same startup -- but this window is not modal, so
            the user could have started working in the main window first.
            Ask before restart() throws that away, same as closing the app
            with unsaved changes already does. """
        if not getattr(self.main.p_session, "changed", False):
            return True
        return self.main.simple_dialog.question(
            "Opening this session discards the current one, which has unsaved changes. Continue?")

    def _on_open(self, widget):
        path = self._selected_path()
        if path is None or not self._confirm_discard_current_session():
            return
        self.window.destroy()
        try:
            self.main.p_session.load_easyhybrid_serialization_file(filename=path)
        except Exception as error:
            self.main.simple_dialog.error(msg="Could not open:\n{}\n\n{}".format(path, error))

    def _on_save_as(self, widget):
        path = self._selected_path()
        if path is None:
            return
        chooser = Gtk.FileChooserDialog(title="Save As...", transient_for=self.window,
                                        action=Gtk.FileChooserAction.SAVE)
        chooser.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        chooser.set_current_name("session.easy")
        chooser.set_do_overwrite_confirmation(True)
        response, dest = chooser.run(), chooser.get_filename()
        chooser.destroy()
        if response != Gtk.ResponseType.OK or not dest:
            return
        if not dest.endswith(".easy"):
            dest += ".easy"
        try:
            shutil.copyfile(path, dest)
        except OSError as error:
            self.main.simple_dialog.error(msg="Could not save:\n{}".format(error))
            return
        self.main.vm_session.vm_config.add_recent_file(dest, "project")
        self.main.simple_dialog.info(msg="Saved to:\n{}".format(dest))

    def _on_delete(self, widget):
        path = self._selected_path()
        if path is None:
            return
        if not self.main.simple_dialog.question("Delete this unsaved session file?\n\n{}".format(path)):
            return
        try:
            os.remove(path)
        except OSError as error:
            self.main.simple_dialog.error(msg="Could not delete:\n{}".format(error))
        self._reload()

    def _on_clear_all(self, widget):
        files = self._list_files()
        if not files:
            return
        if not self.main.simple_dialog.question(
                "Delete all {} unsaved session file(s)? This cannot be undone.".format(len(files))):
            return
        for path, _mtime, _size in files:
            try:
                os.remove(path)
            except OSError:
                pass
        self._reload()

    def _on_destroy(self, widget):
        self.window = None
