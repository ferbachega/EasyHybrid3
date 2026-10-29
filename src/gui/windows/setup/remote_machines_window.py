#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Simulate > Remote Machines... (2026-10-28): the machines reached by ssh
#  where simulations can run (src/remote, docs/REMOTE_JOBS.md). Profiles in
#  ~/.easyhybrid/remotes.json; "Test Connection" checks ssh, the scheduler,
#  Python + pDynamo and computes a real energy there.

import threading

import gi

gi.require_version("Gtk", "3.0")
from gi.repository import GLib, Gtk  # noqa: E402

from remote import profiles as remote_profiles  # noqa: E402

_HELP = ("The host is best an alias of ~/.ssh/config (user, key, port and jump host come from there). "
         "EasyHybrid uses ssh keys only and stores no password. The machine needs pDynamo; "
         "EasyHybrid copies its own code there at each submission.")


class RemoteMachinesWindow:
    def __init__(self, main=None):
        self.main = main
        self.window = None

    # ------------------------------------------------------------------ build
    def open_window(self):
        if self.window is not None:
            self.window.present()
            return
        self.window = Gtk.Window(title="Remote Machines")
        if self.main is not None and getattr(self.main, "window", None) is not None:
            self.window.set_transient_for(self.main.window)
        self.window.set_default_size(740, 900)
        self.window.connect("destroy", self._on_destroy)
        outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8, border_width=10)
        self.window.add(outer)

        top = Gtk.Box(spacing=6)
        top.pack_start(Gtk.Label(label="Machine"), False, False, 0)
        self.combo = Gtk.ComboBoxText()
        self.combo.connect("changed", self._on_select)
        top.pack_start(self.combo, True, True, 0)
        for label, handler in (("New", self._on_new), ("Delete", self._on_delete)):
            button = Gtk.Button(label=label)
            button.connect("clicked", handler)
            top.pack_start(button, False, False, 0)
        outer.pack_start(top, False, False, 0)

        help_label = Gtk.Label(label=_HELP, xalign=0)
        help_label.set_line_wrap(True)
        help_label.get_style_context().add_class("dim-label")
        outer.pack_start(help_label, False, False, 0)

        grid = Gtk.Grid(row_spacing=6, column_spacing=10)
        outer.pack_start(grid, False, False, 0)
        self.fields = {}

        def entry(row, key, label, tooltip, placeholder=""):
            grid.attach(Gtk.Label(label=label, xalign=1), 0, row, 1, 1)
            widget = Gtk.Entry(hexpand=True, placeholder_text=placeholder)
            widget.set_tooltip_text(tooltip)
            grid.attach(widget, 1, row, 3, 1)
            self.fields[key] = widget

        entry(0, "name", "Name", "A short name, shown in Simulate > Run Simulations On.", "bigblue")
        grid.attach(Gtk.Label(label="SSH host", xalign=1), 0, 1, 1, 1)
        hbox = Gtk.Box(spacing=6)
        self.combo_host = Gtk.ComboBoxText.new_with_entry()
        self.combo_host.set_hexpand(True)
        self.combo_host.set_tooltip_text("Alias of ~/.ssh/config (listed), or user@host. The ssh key must already "
                                         "work: 'ssh <host>' in a terminal must log in without a password.")
        from remote.discover import ssh_config_hosts
        for alias in ssh_config_hosts():
            self.combo_host.append_text(alias)
        hbox.pack_start(self.combo_host, True, True, 0)
        self.button_detect = Gtk.Button(label="Detect Settings")
        self.button_detect.set_tooltip_text("Ask the machine (ssh) for its scheduler, partitions, accounts, scratch, "
                                            "the Python that imports pDynamo and the QC programs, and fill the form.")
        self.button_detect.connect("clicked", self._on_detect)
        hbox.pack_start(self.button_detect, False, False, 0)
        button_add_host = Gtk.Button(label="Add SSH Host...")
        button_add_host.set_tooltip_text("A machine known by address, user and port: adds an alias to ~/.ssh/config "
                                         "and shows the commands that authorize your ssh key there.")
        button_add_host.connect("clicked", self._on_add_ssh_host)
        hbox.pack_start(button_add_host, False, False, 0)
        grid.attach(hbox, 1, 1, 3, 1)
        self.fields["host"] = self.combo_host.get_child()
        self.fields["host"].set_placeholder_text("bigblue")
        entry(2, "workdir", "Remote folder", "Code and jobs go here (absolute, or relative to the remote home). "
              "On a cluster, use the scratch file system.", "/scratch/username/easyhybrid")
        entry(3, "python", "Python", "Interpreter that imports pDynamo (after the setup lines).", "python3")

        grid.attach(Gtk.Label(label="Setup lines", xalign=1, yalign=0), 0, 4, 1, 1)
        self.setup_view = Gtk.TextView(monospace=True)
        self.setup_view.set_tooltip_text("Shell lines run before Python on the remote machine: module load, "
                                         "source the pDynamo environment, conda activate...")
        scroller = Gtk.ScrolledWindow(min_content_height=70, hexpand=True)
        scroller.add(self.setup_view)
        grid.attach(scroller, 1, 4, 3, 1)

        grid.attach(Gtk.Label(label="Scheduler", xalign=1), 0, 5, 1, 1)
        self.combo_scheduler = Gtk.ComboBoxText()
        self.combo_scheduler.append("direct", "Direct (background process, workstations)")
        self.combo_scheduler.append("slurm", "SLURM (sbatch)")
        self.combo_scheduler.connect("changed", lambda *_: self._sync_slurm())
        grid.attach(self.combo_scheduler, 1, 5, 1, 1)
        grid.attach(Gtk.Label(label="CPUs per job", xalign=1), 2, 5, 1, 1)
        self.spin_cpus = Gtk.SpinButton.new_with_range(1, 256, 1)
        grid.attach(self.spin_cpus, 3, 5, 1, 1)

        self.slurm_frame = Gtk.Frame(label="SLURM")
        sgrid = Gtk.Grid(row_spacing=6, column_spacing=10, border_width=8)
        self.slurm_frame.add(sgrid)
        self.slurm = {}
        # . partition: the list of the cluster (sinfo), still editable
        sgrid.attach(Gtk.Label(label="Partition", xalign=1), 0, 0, 1, 1)
        pbox = Gtk.Box(spacing=4)
        self.combo_partition = Gtk.ComboBoxText.new_with_entry()
        self.combo_partition.set_hexpand(True)
        self.combo_partition.set_tooltip_text("The partitions of the cluster (the default one first); "
                                              "a name can also be typed.")
        self.combo_partition.connect("changed", self._on_partition_changed)
        pbox.pack_start(self.combo_partition, True, True, 0)
        self.button_partitions = Gtk.Button.new_from_icon_name("view-refresh-symbolic", Gtk.IconSize.BUTTON)
        self.button_partitions.set_tooltip_text("List the partitions of the cluster (sinfo, through ssh)")
        self.button_partitions.connect("clicked", lambda *_: self._load_partitions(force=True))
        pbox.pack_start(self.button_partitions, False, False, 0)
        sgrid.attach(pbox, 1, 0, 1, 1)
        self.slurm["partition"] = self.combo_partition.get_child()
        sgrid.attach(Gtk.Label(label="Account", xalign=1), 2, 0, 1, 1)
        self.combo_account = Gtk.ComboBoxText.new_with_entry()
        self.combo_account.set_hexpand(True)
        self.combo_account.set_tooltip_text("Your allocation, if the cluster asks for one (Detect Settings lists yours).")
        sgrid.attach(self.combo_account, 3, 0, 1, 1)
        self.slurm["account"] = self.combo_account.get_child()
        for i, (key, label, tip) in enumerate((("time", "Time limit", "HH:MM:SS or D-HH:MM:SS"),
                                               ("mem", "Memory", "e.g. 8G (empty: the partition default)")), start=2):
            sgrid.attach(Gtk.Label(label=label, xalign=1), (i % 2) * 2, i // 2, 1, 1)
            widget = Gtk.Entry(hexpand=True)
            widget.set_tooltip_text(tip)
            sgrid.attach(widget, (i % 2) * 2 + 1, i // 2, 1, 1)
            self.slurm[key] = widget
        self.partition_info = Gtk.Label(xalign=0)
        self.partition_info.get_style_context().add_class("dim-label")
        sgrid.attach(self.partition_info, 1, 2, 3, 1)
        sgrid.attach(Gtk.Label(label="Extra #SBATCH", xalign=1), 0, 3, 1, 1)
        self.slurm["extra"] = Gtk.Entry(hexpand=True, placeholder_text="--qos=normal")
        sgrid.attach(self.slurm["extra"], 1, 3, 3, 1)
        self._partitions = {}          # . host -> [partition rows]
        outer.pack_start(self.slurm_frame, False, False, 0)

        entry(6, "fetch_exclude", "Do not fetch", "Space-separated rsync patterns left on the remote machine "
              "when results are fetched, e.g. *.dcd data_collection_traj", "")

        # . pDynamo and the external QC programs on that machine (QC/MM with ORCA, xTB, DFTB+)
        qc_frame = Gtk.Frame(label="pDynamo and QC programs on that machine")
        qgrid = Gtk.Grid(row_spacing=6, column_spacing=10, border_width=8)
        qc_frame.add(qgrid)
        self.qc = {}
        for i, (key, label, tip) in enumerate((
                ("orca", "ORCA", "full path of the orca executable there (parallel ORCA also needs the "
                                 "MPI module it was built with, in the setup lines)"),
                ("xtb", "xTB", "full path of the xtb executable there"),
                ("dftb", "DFTB+", "full path of the dftb+ executable there (the Slater-Koster files are sent "
                                  "with each job)"))):
            qgrid.attach(Gtk.Label(label=label, xalign=1), 0, i, 1, 1)
            widget = Gtk.Entry(hexpand=True, placeholder_text="empty: set by the setup lines, if used")
            widget.set_tooltip_text(tip)
            qgrid.attach(widget, 1, i, 1, 1)
            self.qc[key] = widget
        qgrid.attach(Gtk.Label(label="pDynamo folder", xalign=1), 0, 3, 1, 1)
        self.entry_pdynamo = Gtk.Entry(hexpand=True, placeholder_text="~/programs/pDynamo3")
        self.entry_pdynamo.set_tooltip_text("The pDynamo installation there (its compiled modules suit that machine).")
        qgrid.attach(self.entry_pdynamo, 1, 3, 1, 1)
        self.check_sync = Gtk.CheckButton(label="Use this computer's pDynamo Python files and parameters there")
        self.check_sync.set_tooltip_text("Jobs run on a copy of the pDynamo folder above (in the remote folder) with "
                                         "this computer's .py files and parameters on top: local fixes and DYFF "
                                         "parameters go with the jobs; the installation there is not changed.")
        qgrid.attach(self.check_sync, 1, 4, 1, 1)
        outer.pack_start(qc_frame, False, False, 0)

        buttons = Gtk.Box(spacing=6)
        self.button_test = Gtk.Button(label="Test Connection")
        self.button_test.connect("clicked", self._on_test)
        save = Gtk.Button(label="Save")
        save.connect("clicked", self._on_save)
        buttons.pack_end(save, False, False, 0)
        buttons.pack_end(self.button_test, False, False, 0)
        self.spinner = Gtk.Spinner()
        buttons.pack_start(self.spinner, False, False, 0)
        outer.pack_start(buttons, False, False, 0)

        self.report = Gtk.TextView(editable=False, monospace=True, wrap_mode=Gtk.WrapMode.WORD_CHAR)
        scroller = Gtk.ScrolledWindow(vexpand=True)
        scroller.add(self.report)
        outer.pack_start(scroller, True, True, 0)

        self._reload()
        self.window.show_all()
        self._sync_slurm()

    def _on_destroy(self, *_):
        self.window = None

    # ------------------------------------------------------------------ data
    def _reload(self, select=None):
        self.profiles = remote_profiles.load_all()
        self.combo.remove_all()
        for name in self.profiles:
            self.combo.append(name, name)
        if select and select in self.profiles:
            self.combo.set_active_id(select)
        elif self.profiles:
            self.combo.set_active(0)
        else:
            self._fill(dict(remote_profiles.DEFAULT_PROFILE))

    def _fill(self, p):
        for key in ("name", "host", "workdir", "python"):
            self.fields[key].set_text(str(p.get(key) or ""))
        self.fields["fetch_exclude"].set_text(" ".join(p.get("fetch_exclude") or []))
        self.setup_view.get_buffer().set_text(p.get("setup") or "")
        self.combo_scheduler.set_active_id(p.get("scheduler") or "direct")
        self.spin_cpus.set_value(int(p.get("cpus") or 1))
        for key, widget in self.slurm.items():
            widget.set_text(str((p.get("slurm") or {}).get(key) or ""))
        for key, widget in self.qc.items():
            widget.set_text(str((p.get("qc_programs") or {}).get(key) or ""))
        self.entry_pdynamo.set_text(p.get("pdynamo_home") or "")
        self.check_sync.set_active(bool(p.get("sync_pdynamo")))

    def _read(self):
        buf = self.setup_view.get_buffer()
        return {"name": self.fields["name"].get_text().strip(), "host": self.fields["host"].get_text().strip(),
                "workdir": self.fields["workdir"].get_text().strip(),
                "python": self.fields["python"].get_text().strip() or "python3",
                "setup": buf.get_text(buf.get_start_iter(), buf.get_end_iter(), False).strip(),
                "scheduler": self.combo_scheduler.get_active_id() or "direct",
                "cpus": int(self.spin_cpus.get_value()),
                "slurm": {k: w.get_text().strip() for k, w in self.slurm.items()},
                "fetch_exclude": self.fields["fetch_exclude"].get_text().split(),
                "qc_programs": {k: w.get_text().strip() for k, w in self.qc.items()},
                "pdynamo_home": self.entry_pdynamo.get_text().strip(),
                "sync_pdynamo": self.check_sync.get_active()}

    def _sync_slurm(self):
        slurm = self.combo_scheduler.get_active_id() == "slurm"
        self.slurm_frame.set_sensitive(slurm)
        if slurm:
            self._load_partitions()

    # ------------------------------------------------------------------ partitions
    def _load_partitions(self, force=False):
        """ Fills the partition list from the cluster (sinfo), in a thread;
            once per host unless forced. """
        host = self.fields["host"].get_text().strip()
        if not host or (host in self._partitions and not force):
            self._show_partitions(host)
            return
        self.button_partitions.set_sensitive(False)
        self.partition_info.set_text("Asking {} for its partitions...".format(host))
        profile = self._read()

        def work():
            from remote import schedulers
            from remote import transport as remote_transport
            try:
                rows = schedulers.get("slurm").partitions(remote_transport.make(profile))
                error = None
            except Exception as exc:
                rows, error = None, str(exc)
            GLib.idle_add(self._partitions_loaded, host, rows, error)

        threading.Thread(target=work, daemon=True).start()

    def _partitions_loaded(self, host, rows, error):
        self.button_partitions.set_sensitive(True)
        if rows is None:
            self.partition_info.set_text("Could not list the partitions: {}".format(error)[:160])
            return False
        self._partitions[host] = rows
        self._show_partitions(host)
        return False

    def _show_partitions(self, host):
        rows = self._partitions.get(host)
        if rows is None:
            return
        current = self.slurm["partition"].get_text().strip()
        self.combo_partition.remove_all()
        for row in rows:
            self.combo_partition.append_text(row["name"])
        self.slurm["partition"].set_text(current or next((r["name"] for r in rows if r["default"]), ""))
        self._on_partition_changed(self.combo_partition)

    def _on_partition_changed(self, combo):
        host = self.fields["host"].get_text().strip()
        name = self.slurm["partition"].get_text().strip()
        row = next((r for r in self._partitions.get(host, []) if r["name"] == name), None)
        if row is None:
            self.partition_info.set_text("" if not self._partitions.get(host) else
                                         "'{}' is not a partition of {}.".format(name, host) if name else "")
            return
        text = "{}{} · max {} · {} nodes".format(row["available"], " (default)" if row["default"] else "",
                                                row["time_limit"], row["nodes"])
        if row["cpus_total"]:
            text += " · {} of {} CPUs idle".format(row["cpus_idle"], row["cpus_total"])
        if row["gres"]:
            text += " · {}".format(row["gres"])
        self.partition_info.set_text(text)

    def _say(self, text, append=False):
        buf = self.report.get_buffer()
        if append:
            buf.insert(buf.get_end_iter(), text + "\n")
        else:
            buf.set_text(text + "\n")

    # ------------------------------------------------------------------ actions
    def _on_select(self, combo):
        name = combo.get_active_id()
        if name in getattr(self, "profiles", {}):
            self._fill(self.profiles[name])

    def _on_new(self, _button):
        self.combo.set_active(-1)
        self._fill(dict(remote_profiles.DEFAULT_PROFILE))
        self._say("New machine: choose or type its SSH host, press Detect Settings, review the form, "
                  "then Save and Test Connection.")

    def _on_delete(self, _button):
        name = self.combo.get_active_id()
        if name and remote_profiles.remove(name):
            self._say("Deleted '{}'.".format(name))
            self._reload()
            self._refresh_menu()

    def _on_save(self, _button):
        try:
            saved = remote_profiles.save(self._read())
        except ValueError as error:
            self._say("Not saved: {}.".format(error))
            return
        self._say("Saved '{}'.".format(saved["name"]))
        self._reload(select=saved["name"])
        self._refresh_menu()

    def _refresh_menu(self):
        refresh = getattr(self.main, "refresh_run_on_menu", None)
        if refresh is not None:
            refresh()

    def _on_add_ssh_host(self, _button):
        """ Address + user + port (+ key) -> a Host block in ~/.ssh/config. """
        from remote import sshconfig
        dialog = Gtk.Dialog(title="Add SSH Host", transient_for=self.window, modal=True)
        dialog.add_button("Cancel", Gtk.ResponseType.CANCEL)
        dialog.add_button("Add", Gtk.ResponseType.OK)
        box = dialog.get_content_area()
        box.set_spacing(8)
        box.set_border_width(12)
        grid = Gtk.Grid(row_spacing=6, column_spacing=10)
        box.pack_start(grid, False, False, 0)
        fields = {}
        for row, (key, label, placeholder) in enumerate((("alias", "Alias", "ueg"),
                                                         ("hostname", "Address", "200.137.241.6 or cluster.example.org"),
                                                         ("user", "User", "your login there"),
                                                         ("port", "Port", "22"),
                                                         ("identity", "Key file", "~/.ssh/id_ed25519_<alias>"))):
            grid.attach(Gtk.Label(label=label, xalign=1), 0, row, 1, 1)
            entry = Gtk.Entry(hexpand=True, placeholder_text=placeholder, width_chars=38)
            grid.attach(entry, 1, row, 1, 1)
            fields[key] = entry
        fields["port"].set_text("22")
        preview = Gtk.Label(xalign=0, selectable=True)
        preview.get_style_context().add_class("dim-label")
        box.pack_start(preview, False, False, 0)

        def values():
            alias = fields["alias"].get_text().strip()
            return (alias, fields["hostname"].get_text().strip(), fields["user"].get_text().strip(),
                    fields["port"].get_text().strip() or "22",
                    fields["identity"].get_text().strip() or (sshconfig.default_identity(alias) if alias else ""))

        def update(*_):
            alias, hostname, user, port, identity = values()
            problems = sshconfig.validate(alias, hostname, user, port)
            if not problems and alias in sshconfig.existing_aliases():
                problems = ["'{}' already exists in ~/.ssh/config".format(alias)]
            dialog.set_response_sensitive(Gtk.ResponseType.OK, not problems)
            if problems:
                preview.set_markup("<small>{}</small>".format(GLib.markup_escape_text("; ".join(problems))))
            else:
                text = "Added to ~/.ssh/config (a backup is kept):\n\n" + sshconfig.block(alias, hostname, user, port, identity)
                preview.set_markup("<tt><small>{}</small></tt>".format(GLib.markup_escape_text(text)))

        for entry in fields.values():
            entry.connect("changed", update)
        update()
        dialog.show_all()
        response = dialog.run()
        alias, hostname, user, port, identity = values()
        dialog.destroy()
        if response != Gtk.ResponseType.OK:
            return
        try:
            path, backup = sshconfig.add_host(alias, hostname, user, port, identity)
        except (ValueError, OSError) as error:
            self._say("Not added: {}".format(error))
            return
        self.combo_host.append_text(alias)
        self.fields["host"].set_text(alias)
        if not self.fields["name"].get_text().strip():
            self.fields["name"].set_text(alias)
        commands = sshconfig.setup_commands(alias, hostname, user, port, identity)
        self._say("Added '{}' to {} (backup: {}).\n\nNow, in a terminal (the password of {} is asked only by "
                  "ssh-copy-id, once; EasyHybrid never asks for it):\n\n{}\n\nWhen the last command prints the "
                  "machine's name without asking for a password, press Detect Settings.".format(
                      alias, path, backup or "none", user, "\n".join("  " + c for c in commands)))

    def _on_detect(self, _button):
        host = self.fields["host"].get_text().strip()
        if not host:
            self._say("Choose or type the SSH host first.")
            return
        name = self.fields["name"].get_text().strip() or host.split("@")[-1]
        self.button_detect.set_sensitive(False)
        self.spinner.start()
        self._say("Asking {} for its settings (two ssh calls, up to a minute)...".format(host))

        def work():
            from remote.discover import discover
            try:
                result, error = discover(host, name=name), None
            except Exception as exc:
                result, error = None, str(exc)
            GLib.idle_add(self._detected, host, result, error)

        threading.Thread(target=work, daemon=True).start()

    def _detected(self, host, result, error):
        self.spinner.stop()
        self.button_detect.set_sensitive(True)
        if result is None:
            self._say("Could not ask {}: {}\n\nCheck that 'ssh {}' logs in from a terminal without a password "
                      "(ssh-keygen, then ssh-copy-id {}), and that the VPN is on if the machine needs one.".format(
                          host, error, host, host))
            return False
        profile = result["profile"]
        inv = result["inventory"]
        self._partitions[host] = inv.get("partitions") or []
        self.combo_account.remove_all()
        for account in result["alternatives"].get("accounts", []):
            self.combo_account.append_text(account)
        self._fill(profile)
        self._sync_slurm()
        lines = ["Detected on {} ({}, user {}, {} cores on the login node; system {}):".format(
            host, inv.get("hostname"), inv.get("user"), inv.get("cores"), inv.get("system python"))]
        lines += ["  - " + n for n in result["notes"]] or ["  (nothing to report)"]
        lines.append("\nAlternatives (type them in the form if you prefer one):")
        for key, values in result["alternatives"].items():
            if values:
                lines.append("  {}: {}".format(key, ", ".join(values[:8]) + (" ..." if len(values) > 8 else "")))
        lines.append("\nReview the form, then Save and Test Connection.")
        self._say("\n".join(lines))
        return False

    def _on_test(self, _button):
        profile = self._read()
        problems = remote_profiles.validate(profile)
        if problems:
            self._say("Fix first: {}.".format("; ".join(problems)))
            return
        self.button_test.set_sensitive(False)
        self.spinner.start()
        self._say("Testing {} ({}). The first test copies the EasyHybrid code (~17 MB)...".format(
            profile["name"], profile["host"]))

        def work():
            from remote.jobs import RemoteJobs
            try:
                results = RemoteJobs().test(profile)
            except Exception as error:
                results = [{"check": "test", "ok": False, "detail": "{}: {}".format(type(error).__name__, error)}]
            GLib.idle_add(self._show_test, results)

        threading.Thread(target=work, daemon=True).start()

    def _show_test(self, results):
        self.spinner.stop()
        self.button_test.set_sensitive(True)
        for r in results:
            self._say("{}  {}{}".format("OK  " if r["ok"] else "FAIL", r["check"],
                                        (": " + r["detail"]) if r["detail"] else ""), append=True)
        ok = results and all(r["ok"] for r in results)
        self._say("\nReady: choose it in Simulate > Run Simulations On." if ok else
                  "\nNot ready: see the failed check above.", append=True)
        return False
