#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Edit > AI Assistant (Local Model)... (2026-10-28): a chat with a model
#  served on this computer (Ollama, or an OpenAI-compatible server) that
#  operates this window through the commands of automation/api.py. The model
#  runs in a worker thread (automation/local_agent.py); the commands that
#  touch the interface run in the GTK thread; actions ask for confirmation
#  unless "Run actions without asking" is ticked.

import json
import threading
import time

import gi

gi.require_version("Gtk", "3.0")
from gi.repository import GLib, Gtk, Pango  # noqa: E402

from automation import settings as ai_settings  # noqa: E402


class AIAssistantWindow:
    def __init__(self, main=None):
        self.main = main
        self.window = None
        self.agent = None
        self.busy = False

    # ------------------------------------------------------------------ build
    def open_window(self):
        if self.window is not None:
            self.window.present()
            return
        from automation.local_agent import PROVIDERS
        self.settings = ai_settings.load()
        self.window = Gtk.Window(title="AI Assistant (Local Model)")
        if self.main is not None and getattr(self.main, "window", None) is not None:
            self.window.set_transient_for(self.main.window)
        self.window.set_default_size(820, 640)
        self.window.connect("destroy", self._on_destroy)
        outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, border_width=8)
        self.window.add(outer)

        row = Gtk.Box(spacing=6)
        self.combo_provider = Gtk.ComboBoxText()
        for key, label in PROVIDERS:
            self.combo_provider.append(key, label.split(" (")[0])
        self.combo_provider.set_active_id(self.settings.get("local_provider", "ollama"))
        self.combo_provider.connect("changed", self._on_provider)
        row.pack_start(self.combo_provider, False, False, 0)
        self.entry_url = Gtk.Entry(text=self.settings.get("local_url", "http://localhost:11434"), width_chars=24)
        self.entry_url.set_tooltip_text("Address of the model server on this computer (or your network).")
        row.pack_start(self.entry_url, False, False, 0)
        self.combo_model = Gtk.ComboBoxText.new_with_entry()
        self.combo_model.set_hexpand(True)
        self.combo_model.set_tooltip_text("A model that supports tool calling (Ollama shows 'tools').")
        row.pack_start(self.combo_model, True, True, 0)
        refresh = Gtk.Button.new_from_icon_name("view-refresh-symbolic", Gtk.IconSize.BUTTON)
        refresh.set_tooltip_text("List the models of the server")
        refresh.connect("clicked", lambda *_: self._load_models())
        row.pack_start(refresh, False, False, 0)
        outer.pack_start(row, False, False, 0)

        row2 = Gtk.Box(spacing=12)
        self.check_free = Gtk.CheckButton(label="Run actions without asking")
        self.check_free.set_tooltip_text("Off: before loading files, defining the QC region or starting a "
                                         "simulation, the assistant asks you. Reading the session is always free.")
        self.check_free.set_active(not self.settings.get("local_confirm", True))
        row2.pack_start(self.check_free, False, False, 0)
        row2.pack_start(Gtk.Label(label="Context (tokens)"), False, False, 0)
        self.spin_ctx = Gtk.SpinButton.new_with_range(4096, 131072, 4096)
        self.spin_ctx.set_value(int(self.settings.get("local_num_ctx", 16384)))
        row2.pack_start(self.spin_ctx, False, False, 0)
        note = Gtk.Label(label="Runs on this computer; no data leaves it.", xalign=0)
        note.get_style_context().add_class("dim-label")
        row2.pack_end(note, False, False, 0)
        outer.pack_start(row2, False, False, 0)

        self.view = Gtk.TextView(editable=False, cursor_visible=False, wrap_mode=Gtk.WrapMode.WORD_CHAR,
                                 left_margin=8, right_margin=8, top_margin=6, bottom_margin=6)
        buf = self.view.get_buffer()
        buf.create_tag("user", weight=Pango.Weight.BOLD, foreground="#1f5fa8")
        buf.create_tag("assistant")
        buf.create_tag("tool", foreground="#777777", scale=0.9, family="Monospace")
        buf.create_tag("error", foreground="#b00020")
        buf.create_tag("note", foreground="#777777", style=Pango.Style.ITALIC, scale=0.9)
        buf.create_tag("bold", weight=Pango.Weight.BOLD)
        buf.create_tag("code", family="Monospace")
        scroller = Gtk.ScrolledWindow(vexpand=True)
        scroller.add(self.view)
        outer.pack_start(scroller, True, True, 0)

        bottom = Gtk.Box(spacing=6)
        self.entry = Gtk.Entry(hexpand=True, placeholder_text="Ask, e.g.: what is the barrier of the last scan?")
        self.entry.connect("activate", lambda *_: self._send())
        bottom.pack_start(self.entry, True, True, 0)
        self.button_send = Gtk.Button(label="Send")
        self.button_send.connect("clicked", lambda *_: self._send())
        bottom.pack_start(self.button_send, False, False, 0)
        self.button_stop = Gtk.Button(label="Stop", sensitive=False)
        self.button_stop.connect("clicked", lambda *_: self._stop())
        bottom.pack_start(self.button_stop, False, False, 0)
        clear = Gtk.Button(label="New Conversation")
        clear.connect("clicked", lambda *_: self._clear())
        bottom.pack_start(clear, False, False, 0)
        outer.pack_start(bottom, False, False, 0)
        status = Gtk.Box(spacing=6)
        self.spinner = Gtk.Spinner()
        status.pack_start(self.spinner, False, False, 0)
        self.label_status = Gtk.Label(xalign=0)
        self.label_status.get_style_context().add_class("dim-label")
        status.pack_start(self.label_status, True, True, 0)
        outer.pack_start(status, False, False, 0)

        self.window.show_all()
        self._write("The assistant uses a model running on this computer and the same commands as the "
                    "AI control server (docs/AI_CONTROL.md). Local models answer questions and run "
                    "short tasks well; long protocols are better planned with a larger model.\n", "note")
        self._load_models()

    def _on_destroy(self, *_):
        self._stop()
        self._save_settings()
        self.window = None

    # ------------------------------------------------------------------ helpers
    def _write(self, text, tag="assistant"):
        buf = self.view.get_buffer()
        buf.insert_with_tags_by_name(buf.get_end_iter(), text + ("" if text.endswith("\n") else "\n"), tag)
        mark = buf.create_mark(None, buf.get_end_iter(), False)
        self.view.scroll_mark_onscreen(mark)

    def _write_markdown(self, text):
        """ The model's answer with light Markdown: **bold**, `code`, and
            tables in monospace with aligned columns. """
        import re
        buf = self.view.get_buffer()
        lines = text.splitlines()
        i = 0
        while i < len(lines):
            if lines[i].strip().startswith("|"):
                rows = []
                while i < len(lines) and lines[i].strip().startswith("|"):
                    cells = [c.strip().replace("**", "") for c in lines[i].strip().strip("|").split("|")]
                    if not all(set(c) <= set("-: ") for c in cells):
                        rows.append(cells)
                    i += 1
                widths = [max(len(r[k]) if k < len(r) else 0 for r in rows) for k in range(max(len(r) for r in rows))]
                for n, r in enumerate(rows):
                    line = "  ".join((r[k] if k < len(r) else "").ljust(widths[k]) for k in range(len(widths)))
                    buf.insert_with_tags_by_name(buf.get_end_iter(), line.rstrip() + "\n", "code", *(["bold"] if n == 0 else []))
                continue
            for part in re.split(r"(\*\*[^*]+\*\*|`[^`]+`)", lines[i]):
                if part.startswith("**") and part.endswith("**") and len(part) > 4:
                    buf.insert_with_tags_by_name(buf.get_end_iter(), part[2:-2], "bold")
                elif part.startswith("`") and part.endswith("`") and len(part) > 2:
                    buf.insert_with_tags_by_name(buf.get_end_iter(), part[1:-1], "code")
                elif part:
                    buf.insert_with_tags_by_name(buf.get_end_iter(), part, "assistant")
            buf.insert(buf.get_end_iter(), "\n")
            i += 1
        buf.insert(buf.get_end_iter(), "\n")
        mark = buf.create_mark(None, buf.get_end_iter(), False)
        self.view.scroll_mark_onscreen(mark)

    def _status(self, text):
        self.label_status.set_text(text)

    def _save_settings(self):
        s = ai_settings.load()
        s.update({"local_provider": self.combo_provider.get_active_id() or "ollama",
                  "local_url": self.entry_url.get_text().strip(),
                  "local_model": self.combo_model.get_child().get_text().strip(),
                  "local_num_ctx": int(self.spin_ctx.get_value()),
                  "local_confirm": not self.check_free.get_active()})
        try:
            ai_settings.save(s)
        except OSError:
            pass

    def _on_provider(self, combo):
        from automation.local_agent import DEFAULT_URLS
        self.entry_url.set_text(DEFAULT_URLS.get(combo.get_active_id(), self.entry_url.get_text()))
        self._load_models()

    def _load_models(self):
        provider, url = self.combo_provider.get_active_id() or "ollama", self.entry_url.get_text().strip()
        wanted = self.combo_model.get_child().get_text().strip() or self.settings.get("local_model", "")
        self._status("Listing the models of {}...".format(url))

        def work():
            from automation.local_agent import list_models
            try:
                models, error = list_models(provider, url), None
            except Exception as exc:
                models, error = [], str(exc)
            GLib.idle_add(self._models_loaded, models, error, wanted)

        threading.Thread(target=work, daemon=True).start()

    def _models_loaded(self, models, error, wanted):
        self.combo_model.remove_all()
        for name, note in models:
            self.combo_model.append(name, "{}  ({})".format(name, note) if note else name)
        names = [m for m, _n in models]
        if wanted in names:
            self.combo_model.set_active_id(wanted)
        elif names:
            tools = [m for m, n in models if "tools" in n]
            self.combo_model.set_active_id((tools or names)[0])
        self.combo_model.get_child().set_text(self.combo_model.get_active_id() or wanted)
        self._status(error or "{} model(s) on the server.".format(len(models)))
        if error:
            self._write("Could not list the models: {}\nStart the server (for Ollama: 'ollama serve') and press "
                        "the refresh button.".format(error), "error")
        return False

    # ------------------------------------------------------------------ main-thread calls
    def _in_main(self, fn):
        box, done = {}, threading.Event()

        def run():
            try:
                box["value"] = fn()
            except Exception as error:
                box["error"] = error
            finally:
                done.set()
            return False

        GLib.idle_add(run)
        done.wait()
        if "error" in box:
            raise box["error"]
        return box["value"]

    def _execute(self, name, args):
        """ A command of the API: in the worker thread when it does not touch
            GTK (off_main_thread), else in the GTK thread. """
        from automation.api import AutomationAPI
        if getattr(self, "_api", None) is None:
            settings = dict(ai_settings.load(), summaries=True, compact_json=True)
            self._api = AutomationAPI(self.main, settings=settings)
        fn = getattr(self._api, name)
        if getattr(fn, "off_main_thread", False):
            return fn(**args)
        return self._in_main(lambda: fn(**args))

    def _confirm(self, name, args):
        if self.check_free.get_active():
            return True

        def ask():
            dialog = Gtk.MessageDialog(transient_for=self.window, modal=True, message_type=Gtk.MessageType.QUESTION,
                                       buttons=Gtk.ButtonsType.YES_NO, text="The assistant wants to run: {}".format(name))
            dialog.format_secondary_text(json.dumps(args, indent=1, default=str)[:1500] + "\n\nAllow it?")
            answer = dialog.run() == Gtk.ResponseType.YES
            dialog.destroy()
            return answer

        allowed = self._in_main(ask)
        GLib.idle_add(self._write, "  {} {}".format("allowed:" if allowed else "refused:", name), "note")
        return allowed

    # ------------------------------------------------------------------ conversation
    def _event(self, kind, data):
        def show():
            if kind == "model":
                self._status("model answered in {:.0f} s{}".format(
                    data["seconds"], (": calls " + ", ".join(data["tools"])) if data["tools"] else ""))
            elif kind == "thought" and data["text"].strip():
                self._write(data["text"].strip(), "note")
            elif kind == "tool_call":
                self._write("  → {}({})".format(data["name"], json.dumps(data["arguments"], default=str)[:300]), "tool")
            elif kind == "tool_result":
                self._write("  ← {}{}".format("" if data["ok"] else "error: ", data["text"][:300]), "tool")
            elif kind == "answer":
                self._write_markdown(data["text"].strip())
            return False
        GLib.idle_add(show)

    def _send(self):
        text = self.entry.get_text().strip()
        if not text or self.busy:
            return
        model = self.combo_model.get_child().get_text().strip()
        if not model:
            self._write("Choose a model first.", "error")
            return
        self._save_settings()
        from automation.local_agent import LocalAgent
        provider, url = self.combo_provider.get_active_id() or "ollama", self.entry_url.get_text().strip()
        if self.agent is None or (self.agent.model, self.agent.url, self.agent.provider) != (model, url.rstrip("/"), provider):
            self.agent = LocalAgent(provider=provider, url=url, model=model, execute=self._execute,
                                    confirm=self._confirm, on_event=self._event, num_ctx=int(self.spin_ctx.get_value()))
        self.agent.num_ctx = int(self.spin_ctx.get_value())
        self.entry.set_text("")
        self._write("You: " + text, "user")
        self.busy = True
        self.button_send.set_sensitive(False)
        self.button_stop.set_sensitive(True)
        self.spinner.start()
        self._status("The model is working (the first answer also loads it)...")
        started = time.time()

        def work():
            try:
                self.agent.ask(text)
                error = None
            except Exception as exc:
                error = str(exc)
            GLib.idle_add(self._done, error, time.time() - started)

        threading.Thread(target=work, daemon=True).start()

    def _done(self, error, seconds):
        self.busy = False
        self.spinner.stop()
        self.button_send.set_sensitive(True)
        self.button_stop.set_sensitive(False)
        if error:
            self._write("Error: " + error, "error")
        self._status("Done in {:.0f} s.".format(seconds) if not error else "Stopped by an error.")
        return False

    def _stop(self):
        if self.agent is not None:
            self.agent.stopped = True
        self._status("Stopping after the current step...")

    def _clear(self):
        if self.busy:
            return
        self.agent = None
        self.view.get_buffer().set_text("")
        self._status("New conversation.")
