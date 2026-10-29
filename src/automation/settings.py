#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Settings of the AI control layer (2026-10-28), chosen in Edit > AI Control
#  Server and saved in ~/.easyhybrid/ai_settings.json (not in .config.json:
#  the MCP server reads them without EasyHybrid running, and several
#  EasyHybrid windows may be open). Standard library only.
#
#    mode              "economical" or "detailed": the presets below
#    summaries         long tables and matrices left out of the replies (the
#                      files they come from are named; detail="full" brings them)
#    auto_wait         optimize/scan/md/neb/umbrella wait for the job and
#                      return its result in the same call
#    wait_limit        seconds a call waits before returning "still running"
#    screenshot_width  screenshots scaled down to this width (0: full size)
#    short_tools       one-line tool descriptions in the MCP tool list (the
#                      "help" tool gives the full text); new Claude Code sessions
#    compact_json      replies without indentation
#    dialogs           information/"Finished!" dialogs still shown
#    operator_model    model of the Claude Code subagent that runs the commands:
#                      "session" (no subagent), "haiku", "sonnet" or "opus"

import json
import os

CONTROL_DIR = os.path.join(os.path.expanduser("~"), ".easyhybrid")
SETTINGS_FILE = os.path.join(CONTROL_DIR, "ai_settings.json")
AGENT_FILE = os.path.join(os.path.expanduser("~"), ".claude", "agents", "easyhybrid-operator.md")
AGENT_MARK = "<!-- written by EasyHybrid (Edit > AI Control Server); changes here are overwritten -->"

PRESETS = {
    "economical": {"summaries": True, "auto_wait": True, "wait_limit": 600, "screenshot_width": 1000,
                   "short_tools": True, "compact_json": True},
    "detailed":   {"summaries": False, "auto_wait": False, "wait_limit": 600, "screenshot_width": 0,
                   "short_tools": False, "compact_json": False},
}
MODELS = (("session", "Same as the Claude Code session (no subagent)"),
          ("haiku", "Claude Haiku: cheapest and fastest, for routine runs"),
          ("sonnet", "Claude Sonnet: balanced"),
          ("opus", "Claude Opus: most capable, most expensive"))

DEFAULTS = dict(PRESETS["economical"], mode="economical", dialogs=False, operator_model="session",
                # . Edit > AI Assistant (local model): automation/local_agent.py
                local_provider="ollama", local_url="http://localhost:11434", local_model="",
                local_num_ctx=16384, local_confirm=True)


def load(path=None):
    """ The saved settings over the defaults (unknown keys ignored). """
    out = dict(DEFAULTS)
    try:
        with open(path or SETTINGS_FILE) as f:
            saved = json.load(f)
        out.update({k: v for k, v in saved.items() if k in DEFAULTS})
    except (OSError, ValueError):
        pass
    return out


def save(settings, path=None):
    path = path or SETTINGS_FILE
    os.makedirs(os.path.dirname(path), mode=0o700, exist_ok=True)
    clean = {k: settings[k] for k in DEFAULTS if k in settings}
    with open(path, "w") as f:
        json.dump(clean, f, indent=1)
    return clean


def apply_preset(settings, mode):
    """ settings with the values of a preset ("economical"/"detailed"). """
    out = dict(settings)
    out.update(PRESETS[mode])
    out["mode"] = mode
    return out


AGENT_TEXT = """---
name: easyhybrid-operator
description: Runs EasyHybrid through its MCP tools (mcp__easyhybrid__*): loads systems, starts and waits for optimizations, scans, MD, NEB and umbrella sampling, and returns a short summary with the numbers. Use it for the execution steps of EasyHybrid work; decide the protocol and interpret the results in the main session.
model: {model}
---
{mark}

You operate EasyHybrid (QC/MM and molecular simulations, pDynamo3) through the `easyhybrid` MCP tools.

- Start with `state` (systems = sys ids, objects = obj ids). Atom indices are 0-based.
- Simulations: optimize, scan, md, neb, umbrella. With auto-wait on they return the finished
  result; otherwise call `wait`, then `import_result` (not for umbrella), then `pes`.
- Umbrella sampling: after it finishes, `us_check`; report its verdict, and do not report a PMF
  as reliable when the verdict is "fail". Then `wham`.
- A reply with "pending" means a dialog is open: `dialogs`, then `dialog_respond`, then `result`.
- Use the `help` tool for the arguments of a command. Prefer numbers to screenshots.
- Energies are in kJ/mol, distances in A. Write "QC/MM", never "QM/MM".

Reply with a short summary: what was run, the ids of new objects and jobs, the key numbers,
and any warning from the tools. Do not paste whole tool outputs.
"""


def write_operator_agent(model, path=None):
    """ Writes (model haiku/sonnet/opus) or removes ("session") the Claude Code
        subagent. Only a file EasyHybrid wrote itself is overwritten or removed.
        Returns the path written, None when removed or left alone. """
    path = path or AGENT_FILE
    ours = True
    if os.path.exists(path):
        try:
            with open(path) as f:
                ours = AGENT_MARK in f.read()
        except OSError:
            ours = False
    if not ours:
        return None
    if model == "session":
        if os.path.exists(path):
            os.unlink(path)
        return None
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(AGENT_TEXT.format(model=model, mark=AGENT_MARK))
    return path
