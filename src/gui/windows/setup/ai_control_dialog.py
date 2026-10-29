#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Edit > AI Control Server (2026-10-28): options of the control server for
#  scripts and AI assistants (src/automation, docs/AI_CONTROL.md) -- token
#  use, the model that operates EasyHybrid in Claude Code, dialogs.
#  Saved in ~/.easyhybrid/ai_settings.json (automation/settings.py).

import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk  # noqa: E402

from automation import settings as ai_settings  # noqa: E402

_MODE_TEXT = {
    "economical": "Short replies, results in the same call, small screenshots: several times fewer tokens.",
    "detailed": "Full data tables and full-size screenshots in every reply: more tokens, nothing left out.",
}


class AIControlDialog:
    """ Modal dialog; run() returns the new settings dict, or None on Cancel. """

    def __init__(self, parent=None, settings=None, starting=True):
        self.settings = dict(settings or ai_settings.load())
        self._updating = False
        self.dialog = Gtk.Dialog(title="AI Control Server", transient_for=parent, modal=True)
        self.dialog.add_button("Cancel", Gtk.ResponseType.CANCEL)
        self.dialog.add_button("Start" if starting else "Apply", Gtk.ResponseType.OK)
        self.dialog.set_default_response(Gtk.ResponseType.OK)
        self.dialog.set_default_size(560, -1)

        box = self.dialog.get_content_area()
        box.set_spacing(10)
        box.set_border_width(12)

        intro = Gtk.Label(xalign=0)
        intro.set_line_wrap(True)
        intro.set_max_width_chars(70)
        intro.set_markup("Lets scripts and AI assistants (Claude Code through the MCP server, or "
                         "<tt>ehctl</tt>) drive this window through a local socket that only your user "
                         "can open. See <tt>docs/AI_CONTROL.md</tt>.")
        box.pack_start(intro, False, False, 0)

        # ------------------------------------------------------------ tokens
        frame = Gtk.Frame(label="Token use")
        grid = Gtk.Grid(row_spacing=6, column_spacing=12, border_width=10)
        frame.add(grid)
        box.pack_start(frame, False, False, 0)

        self.radio_economical = Gtk.RadioButton.new_with_label(None, "Economical (recommended)")
        self.radio_detailed = Gtk.RadioButton.new_with_label_from_widget(self.radio_economical, "Detailed")
        self.mode_label = Gtk.Label(xalign=0)
        self.mode_label.set_line_wrap(True)
        self.mode_label.set_max_width_chars(60)
        self.mode_label.get_style_context().add_class("dim-label")
        grid.attach(self.radio_economical, 0, 0, 1, 1)
        grid.attach(self.radio_detailed, 1, 0, 1, 1)
        grid.attach(self.mode_label, 0, 1, 2, 1)

        self.check_summaries = Gtk.CheckButton(label="Summaries instead of full data tables")
        self.check_summaries.set_tooltip_text("Scan surfaces, PMFs and window tables are left out of the replies "
                                              "(the assistant can still ask for them with detail=\"full\").")
        self.check_wait = Gtk.CheckButton(label="Simulations return their result in the same call")
        self.check_wait.set_tooltip_text("optimize, scan, md, neb and umbrella wait for the job, import it and "
                                         "summarize it, instead of one call for each step.")
        self.check_short = Gtk.CheckButton(label="Short tool descriptions (new Claude Code sessions)")
        self.check_short.set_tooltip_text("One line per tool in the list the assistant loads at the start of a "
                                          "session; the full text is given by the help tool when needed.")
        self.check_compact = Gtk.CheckButton(label="Compact JSON replies")
        width_box = Gtk.Box(spacing=6)
        width_box.pack_start(Gtk.Label(label="Screenshot width (px, 0 = full size)"), False, False, 0)
        self.spin_width = Gtk.SpinButton.new_with_range(0, 3840, 100)
        width_box.pack_start(self.spin_width, False, False, 0)
        wait_box = Gtk.Box(spacing=6)
        wait_box.pack_start(Gtk.Label(label="Wait at most (s) before returning \"still running\""), False, False, 0)
        self.spin_wait = Gtk.SpinButton.new_with_range(30, 7200, 60)
        wait_box.pack_start(self.spin_wait, False, False, 0)
        for row, widget in enumerate((self.check_summaries, self.check_wait, wait_box, self.check_short,
                                      self.check_compact, width_box), start=2):
            grid.attach(widget, 0, row, 2, 1)

        # ------------------------------------------------------------- model
        frame = Gtk.Frame(label="Model that operates EasyHybrid")
        mbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, border_width=10)
        frame.add(mbox)
        box.pack_start(frame, False, False, 0)
        self.combo_model = Gtk.ComboBoxText()
        for key, label in ai_settings.MODELS:
            self.combo_model.append(key, label)
        mbox.pack_start(self.combo_model, False, False, 0)
        note = Gtk.Label(xalign=0)
        note.set_line_wrap(True)
        note.set_max_width_chars(70)
        note.get_style_context().add_class("dim-label")
        note.set_markup("With a model chosen, EasyHybrid writes the Claude Code subagent "
                        "<tt>~/.claude/agents/easyhybrid-operator.md</tt>: it runs the commands with that model "
                        "while the main session plans and interprets. The model of the main session is chosen "
                        "in Claude Code itself (<tt>/model</tt>). Applies to new Claude Code sessions.")
        mbox.pack_start(note, False, False, 0)

        # ----------------------------------------------------------- dialogs
        self.check_dialogs = Gtk.CheckButton(label="Show information and \"Finished!\" dialogs while controlled")
        self.check_dialogs.set_tooltip_text("Off: their messages go to the assistant with the reply. Yes/No "
                                            "questions are always shown.")
        box.pack_start(self.check_dialogs, False, False, 0)

        self._load()
        self.radio_economical.connect("toggled", self._on_mode)
        for widget in (self.check_summaries, self.check_wait, self.check_short, self.check_compact):
            widget.connect("toggled", self._on_custom)
        self.spin_width.connect("value-changed", self._on_custom)
        self.spin_wait.connect("value-changed", self._on_custom)
        self.dialog.show_all()

    # ------------------------------------------------------------------
    def _load(self):
        self._updating = True
        s = self.settings
        (self.radio_economical if s.get("mode") == "economical" else self.radio_detailed).set_active(True)
        self.check_summaries.set_active(bool(s.get("summaries")))
        self.check_wait.set_active(bool(s.get("auto_wait")))
        self.check_short.set_active(bool(s.get("short_tools")))
        self.check_compact.set_active(bool(s.get("compact_json")))
        self.spin_width.set_value(int(s.get("screenshot_width") or 0))
        self.spin_wait.set_value(int(s.get("wait_limit") or 600))
        self.combo_model.set_active_id(s.get("operator_model", "session"))
        self.check_dialogs.set_active(bool(s.get("dialogs")))
        self._show_mode()
        self._updating = False

    def _show_mode(self):
        mode = self.settings.get("mode")
        self.mode_label.set_text(_MODE_TEXT.get(mode, "Custom: the options below."))

    def _on_mode(self, _widget):
        if self._updating:
            return
        mode = "economical" if self.radio_economical.get_active() else "detailed"
        self.settings = ai_settings.apply_preset(self._read(), mode)
        self._load()

    def _on_custom(self, _widget):
        if self._updating:
            return
        current = self._read()
        for mode, preset in ai_settings.PRESETS.items():
            if all(current.get(k) == v for k, v in preset.items()):
                current["mode"] = mode
                break
        else:
            current["mode"] = "custom"
        self.settings = current
        self._show_mode()

    def _read(self):
        out = dict(self.settings)
        out.update({"summaries": self.check_summaries.get_active(), "auto_wait": self.check_wait.get_active(),
                    "short_tools": self.check_short.get_active(), "compact_json": self.check_compact.get_active(),
                    "screenshot_width": int(self.spin_width.get_value()), "wait_limit": int(self.spin_wait.get_value()),
                    "operator_model": self.combo_model.get_active_id() or "session",
                    "dialogs": self.check_dialogs.get_active()})
        return out

    def run(self):
        response = self.dialog.run()
        result = self._on_ok() if response == Gtk.ResponseType.OK else None
        self.dialog.destroy()
        return result

    def _on_ok(self):
        settings = self._read()
        settings["mode"] = self.settings.get("mode", "custom")
        return settings


def save_and_apply(settings):
    """ Saves the settings and writes/removes the Claude Code subagent; returns
        a status line for the main statusbar. """
    ai_settings.save(settings)
    model = settings.get("operator_model", "session")
    try:
        path = ai_settings.write_operator_agent(model)
    except OSError as error:
        return "AI settings saved; could not write the Claude Code subagent: {}".format(error)
    if model == "session":
        return "AI settings saved ({} token use).".format(settings.get("mode"))
    if path is None:
        return ("AI settings saved; ~/.claude/agents/easyhybrid-operator.md was not written by EasyHybrid "
                "and was left unchanged.")
    return "AI settings saved ({} token use, operator model: {}).".format(settings.get("mode"), model)
