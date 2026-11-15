"""
Testes da camada de automacao (src/automation): protocolo do servidor de
controle, cliente, esquema dos comandos e servidor MCP. Sem interface: o
"main" e' um objeto falso e um GLib.MainLoop roda numa thread no lugar do
loop do GTK.
"""
import io
import json
import os
import threading

import pytest

from automation import events
from automation.api import _arg_docs, command_schemas, parse_indices

gi = pytest.importorskip("gi")
from gi.repository import GLib  # noqa: E402


class _FakeMain:
    version = "test"


@pytest.fixture
def server(tmp_path, monkeypatch):
    from automation import server as server_module
    monkeypatch.setattr(server_module, "INFO_FILE", str(tmp_path / "control.json"))
    loop = GLib.MainLoop()
    thread = threading.Thread(target=loop.run, daemon=True)
    thread.start()
    srv = server_module.ControlServer(_FakeMain(), path=str(tmp_path / "eh.sock"), exit_handlers=False)
    srv.start()
    yield srv
    srv.stop()
    loop.quit()
    events.set_active(False)


def test_parse_indices():
    assert parse_indices("0-3,7") == [0, 1, 2, 3, 7]
    assert parse_indices([5, "2"]) == [5, 2]
    assert parse_indices(4) == [4]
    assert parse_indices(None) is None


def test_arg_docs_multiple_names():
    doc = """Summary.

    Args:
      obj, frame: starting coordinates
      name: trajectory name
        continued here
    """
    docs = _arg_docs(doc)
    assert docs["obj"] == docs["frame"] == "starting coordinates"
    assert docs["name"] == "trajectory name continued here"


def test_command_schemas_complete():
    schemas = {s["name"]: s for s in command_schemas()}
    for name in ("state", "atoms", "measure", "set_qc", "optimize", "scan", "md", "neb",
                 "jobs", "import_result", "pes", "screenshot", "dialogs", "dialog_respond", "python",
                 "umbrella", "us_check", "wham", "us_browser"):
        assert name in schemas
    for s in schemas.values():
        assert s["summary"], s["name"]
        for arg, prop in s["inputSchema"]["properties"].items():
            assert "type" in prop, (s["name"], arg)
    assert schemas["scan"]["inputSchema"]["required"] == ["rc1"]


def test_events_since():
    first = events.record("info", "a")
    events.record("error", "b")
    kinds = [e["kind"] for e in events.since(first - 1)]
    assert kinds[-2:] == ["info", "error"]
    assert [e["text"] for e in events.since(first, ["error"])][-1] == "b"


def test_round_trip(server):
    from automation.client import EasyHybridClient, EasyHybridError
    with EasyHybridClient(server.path, timeout=10) as client:
        assert client.call("ping")["version"] == "test"
        reply = client.request("no_such_command")
        assert reply["ok"] is False and "Unknown command" in reply["error"]
        reply = client.request("ping", {"unexpected": 1})
        assert reply["ok"] is False and "arguments" in reply["error"]
        assert any(s["name"] == "scan" for s in client.call("describe"))
        seq = events.record("finished", "job done", system=0, job=1)
        out = client.call("events", since=seq - 1)
        assert out["events"][0]["text"] == "job done"
        with pytest.raises(EasyHybridError):
            client.call("_system")
    assert oct(os.stat(server.path).st_mode & 0o777) == "0o600"
    assert events.active() is True


def test_events_returned_with_reply(server):
    from automation.client import EasyHybridClient
    with EasyHybridClient(server.path, timeout=10) as client:
        reply = client.request("python", {"code": "from automation import events\nevents.record('info', 'hello')\nresult = 41 + 1"})
    assert reply["ok"] and reply["result"]["result"] == 42
    assert [e["text"] for e in reply["events"]] == ["hello"]


def test_mcp_server(server, monkeypatch):
    from automation.mcp_server import MCPServer
    monkeypatch.setenv("EASYHYBRID_SOCKET", server.path)
    mcp = MCPServer()
    init = mcp.handle({"jsonrpc": "2.0", "id": 1, "method": "initialize",
                       "params": {"protocolVersion": "2025-06-18"}})
    assert init["result"]["serverInfo"]["name"] == "easyhybrid"
    assert mcp.handle({"jsonrpc": "2.0", "method": "notifications/initialized"}) is None
    tools = mcp.handle({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})["result"]["tools"]
    names = {t["name"] for t in tools}
    assert {"state", "launch", "wait", "result"} <= names
    call = mcp.handle({"jsonrpc": "2.0", "id": 3, "method": "tools/call",
                       "params": {"name": "ping", "arguments": {}}})["result"]
    assert call["isError"] is False
    assert json.loads(call["content"][0]["text"])["result"]["version"] == "test"
    bad = mcp.handle({"jsonrpc": "2.0", "id": 4, "method": "nope"})
    assert bad["error"]["code"] == -32601
    # . stdio framing: one JSON per line
    out = io.StringIO()
    mcp.serve(io.StringIO('{"jsonrpc":"2.0","id":9,"method":"ping"}\n'), out)
    assert json.loads(out.getvalue()) == {"jsonrpc": "2.0", "id": 9, "result": {}}


def test_remove_stale_sockets(tmp_path):
    from automation.server import ControlServer
    dead = tmp_path / "control-999999999.sock"
    alive = tmp_path / "control-{}.sock".format(os.getpid())
    other = tmp_path / "notes.txt"
    for f in (dead, alive, other):
        f.write_text("")
    ControlServer._remove_stale_sockets(str(tmp_path))
    assert not dead.exists() and alive.exists() and other.exists()


def test_terminal_positional_words():
    pytest.importorskip("OpenGL")
    from gui.windows.setup.easyhybrid_terminal import Command

    class Probe(Command):
        def cmd_echo(self, rep="lines", obj=None, **_):
            return "{} {}".format(rep, obj)

    cmd = Probe(None)
    assert cmd.run("echo sticks") == "sticks None"
    assert cmd.run("echo sticks 3") == "sticks 3"
    assert cmd.run("echo obj=2 spheres") == "spheres 2"
    assert cmd.run("echo a b c").startswith("Usage: echo")


def test_profile_prominence():
    import numpy as np
    from automation.api import _profile_1d
    x = np.linspace(0.0, 4.0, 81)
    F = 10.0 * np.exp(-((x - 2.0) / 0.4) ** 2) + 1.5 * x             # . barrier at 2, uphill end
    F[10] += 0.3                                                        # . a noise ripple
    out = _profile_1d(x, F)
    assert [round(m["rc"], 1) for m in out["local_maxima"]] == [2.0]
    assert out["highest"]["rc"] == pytest.approx(2.0, abs=0.06)
    assert len(out["local_minima"]) == 1 and out["local_minima"][0]["rc"] > 2.0


class _Model:
    def __init__(self, target, k):
        self.lowEquilibriumValue, self.lowForceConstant = target, k


def _synthetic_us(folder, shift_window=None, k=500.0, n=4000):
    """ Flat free energy: each window samples its own restraint, N(target, sqrt(kT/2k)). """
    import pickle
    import numpy as np
    rng = np.random.default_rng(7)
    sigma = (0.0083144626 * 300.0 / (2.0 * k)) ** 0.5
    targets = np.arange(1.0, 2.0, 0.05)
    (folder / "output.log").write_text("no parameters\n")
    for i, t in enumerate(targets):
        w = folder / "data_collection" / "window{}.ptRes".format(i)
        w.mkdir(parents=True)
        with open(w / "header.pkl", "wb") as f:
            pickle.dump({"labels": ["RC1"], "restraintEnergyModels": [_Model(float(t), k)]}, f)
        mean = t + (5 * sigma if i == shift_window else 0.0)
        with open(w / "block0.pkl", "wb") as f:
            pickle.dump(list(rng.normal(mean, sigma, n)), f)
    return sigma


def test_us_check_synthetic(tmp_path):
    from automation.api import AutomationAPI
    good, bad = tmp_path / "good", tmp_path / "bad"
    good.mkdir(), bad.mkdir()
    _synthetic_us(good)
    _synthetic_us(bad, shift_window=7)
    api = AutomationAPI(None)
    ok = api.us_check(folder=str(good))
    assert ok["verdict"] == "OK" and ok["windows_found"] == 20 and not ok["deviating_windows"]
    assert ok["min_overlap"] > 0.03 and ok["halves_difference_within_targets"] < 4.0
    fail = api.us_check(folder=str(bad), halves=False)
    assert fail["verdict"] == "fail" and fail["deviating_windows"] == [[7]]


def test_off_main_thread_commands_run_in_a_worker(server):
    from automation.api import off_main_thread
    from automation.client import EasyHybridClient

    @off_main_thread
    def where():
        return threading.current_thread().name

    server.api.where = where
    server.api.here = lambda: threading.current_thread().name
    with EasyHybridClient(server.path, timeout=10) as client:
        assert client.call("where") == "easyhybrid-command"
        assert client.call("here") != "easyhybrid-command"


def test_edit_menu_toggle(tmp_path, monkeypatch):
    """ Edit > AI Control Server starts/stops the server, and follows a server
        started or stopped elsewhere (--control, terminal 'control'). """
    pytest.importorskip("OpenGL")
    gi.require_version("Gtk", "3.0")
    from gi.repository import Gtk
    from automation import server as server_module
    from gui.main.main_window import MainWindow
    monkeypatch.setattr(server_module, "INFO_FILE", str(tmp_path / "control.json"))
    monkeypatch.setattr(server_module.ControlServer, "_install_exit_handlers", lambda self: None)
    monkeypatch.setenv("EASYHYBRID_SOCKET", str(tmp_path / "m.sock"))
    from automation import settings as ai_settings
    from gui.windows.setup import ai_control_dialog
    monkeypatch.setattr(ai_settings, "SETTINGS_FILE", str(tmp_path / "ai_settings.json"))
    monkeypatch.setattr(ai_settings, "AGENT_FILE", str(tmp_path / "agents" / "easyhybrid-operator.md"))
    chosen = dict(ai_settings.DEFAULTS, operator_model="haiku")
    monkeypatch.setattr(ai_control_dialog.AIControlDialog, "__init__", lambda self, **kw: None)
    monkeypatch.setattr(ai_control_dialog.AIControlDialog, "run", lambda self: dict(chosen))

    class Main:
        version = "test"
        window = None
        on_menuitem_control_server_toggled = MainWindow.on_menuitem_control_server_toggled
        sync_control_menu = MainWindow.sync_control_menu

        def __init__(self):
            self.item = Gtk.CheckMenuItem(label="AI Control Server")
            self.item.connect("toggled", self.on_menuitem_control_server_toggled)
            self.statusbar_main = Gtk.Statusbar()
            self.builder = type("B", (), {"get_object": lambda _b, name: self.item})()

    main = Main()
    main.item.set_active(True)                                    # . menu on -> dialog -> server on
    assert main.control_server.running and os.path.exists(str(tmp_path / "m.sock"))
    assert main.control_server.settings["operator_model"] == "haiku"
    assert "model: haiku" in (tmp_path / "agents" / "easyhybrid-operator.md").read_text()
    assert json.loads((tmp_path / "ai_settings.json").read_text())["operator_model"] == "haiku"
    main.item.set_active(False)                                   # . menu off -> server off
    assert not main.control_server.running and not os.path.exists(str(tmp_path / "m.sock"))
    monkeypatch.setattr(ai_control_dialog.AIControlDialog, "run", lambda self: None)
    main.item.set_active(True)                                    # . dialog cancelled -> stays off
    assert main.item.get_active() is False and not main.control_server.running
    srv = server_module.ControlServer(main, exit_handlers=False)  # . started elsewhere
    main.control_server = srv
    srv.start()
    assert main.item.get_active() is True
    srv.stop()
    assert main.item.get_active() is False
    events.set_active(False)



def test_settings_presets_and_agent(tmp_path):
    from automation import settings as ai_settings
    path = str(tmp_path / "ai.json")
    assert ai_settings.load(path) == ai_settings.DEFAULTS
    detailed = ai_settings.apply_preset(ai_settings.DEFAULTS, "detailed")
    assert detailed["summaries"] is False and detailed["mode"] == "detailed"
    ai_settings.save(dict(detailed, unknown=1), path)
    assert ai_settings.load(path)["mode"] == "detailed" and "unknown" not in ai_settings.load(path)
    agent = str(tmp_path / "agents" / "op.md")
    assert ai_settings.write_operator_agent("sonnet", agent) == agent
    text = open(agent).read()
    assert "model: sonnet" in text and ai_settings.AGENT_MARK in text
    assert ai_settings.write_operator_agent("session", agent) is None and not os.path.exists(agent)
    with open(agent, "w") as f:                       # . a subagent the user wrote: never touched
        f.write("---\nname: easyhybrid-operator\nmodel: opus\n---\nmine\n")
    assert ai_settings.write_operator_agent("haiku", agent) is None
    assert open(agent).read().endswith("mine\n")


def test_auto_wait_and_compact_replies(server):
    from automation.client import EasyHybridClient
    calls = []
    server.api.optimize = lambda **kw: calls.append("optimize") or {"system": 0, "job": 3, "type": "Geometry_Optimization"}
    server.api.job_status = lambda sys=None, job=None: {"status": "Finished", "warning": None}
    server.api.import_result = lambda sys=None, job=None: calls.append("import") or {"object": 7, "note": None}
    server.settings.update(auto_wait=True, compact_json=True)
    with EasyHybridClient(server.path, timeout=10) as client:
        reply = client.request("optimize", {})
        assert reply["format"] == "compact" and "events" not in reply          # . empty lists dropped
        assert reply["result"]["status"] == "Finished" and reply["result"]["imported"] == {"object": 7}
        assert calls == ["optimize", "import"]
        reply = client.request("optimize", {"wait": False})
        assert reply["result"] == {"system": 0, "job": 3, "type": "Geometry_Optimization"}


def test_summary_detail():
    import numpy as np
    from automation.api import AutomationAPI, _drop_arrays
    api = AutomationAPI(None, settings={"summaries": True})
    assert api._summary(None) and not api._summary("full")
    out = _drop_arrays({"RC1": [[1, 2]], "Z": list(range(50)), "path": {"cells": [[0, 0], [1, 1]]}},
                       ("RC1", "Z"), keep_short=True)
    assert "RC1" not in out and "Z" not in out and out["path"]["cells"] == 2 and "omitted" in out
    out = _drop_arrays({"Z": list(np.arange(5.0))}, ("Z",), keep_short=True)
    assert "Z" in out and "omitted" not in out


def test_mcp_short_tools_and_instructions():
    from automation.mcp_server import MCPServer, instructions, tool_list
    from automation import settings as ai_settings
    short = tool_list(ai_settings.DEFAULTS)
    long = tool_list(ai_settings.apply_preset(ai_settings.DEFAULTS, "detailed"))
    assert len(json.dumps(short)) < 0.75 * len(json.dumps(long))
    assert "help" in {t["name"] for t in short}
    text = instructions(dict(ai_settings.DEFAULTS, operator_model="haiku"))
    assert "easyhybrid-operator" in text and "Haiku" in text
    assert "easyhybrid-operator" not in instructions(ai_settings.DEFAULTS)
    reply = MCPServer().handle({"jsonrpc": "2.0", "id": 1, "method": "tools/call",
                                "params": {"name": "help", "arguments": {"tool": "scan"}}})["result"]
    assert not reply["isError"] and "multiple_distance" in reply["content"][0]["text"]


def test_ai_control_dialog_presets():
    pytest.importorskip("OpenGL")
    from gui.windows.setup.ai_control_dialog import AIControlDialog
    from automation import settings as ai_settings
    dialog = AIControlDialog(settings=ai_settings.DEFAULTS)
    assert dialog.radio_economical.get_active() and dialog.check_summaries.get_active()
    dialog.radio_detailed.set_active(True)
    assert dialog.settings["mode"] == "detailed" and not dialog.check_summaries.get_active()
    assert dialog.spin_width.get_value() == 0
    dialog.check_wait.set_active(True)
    assert dialog.settings["mode"] == "custom"
    dialog.combo_model.set_active_id("sonnet")
    assert dialog._on_ok()["operator_model"] == "sonnet"
    dialog.dialog.destroy()


def test_neighbour_overlaps_from_files(tmp_path):
    """ The overlap the WHAM tool now checks (pdynamo.us_data.neighbour_overlaps)
        matches the Browser's, and k = 4000 with 0.05 A spacing gives gaps. """
    import glob
    from pdynamo.us_data import USRun, neighbour_overlaps
    good, stiff = tmp_path / "good", tmp_path / "stiff"
    good.mkdir(), stiff.mkdir()
    _synthetic_us(good)                          # . k = 500: sigma 0.050, spacing 0.05
    _synthetic_us(stiff, k=40000.0)              # . sigma 0.0056: windows 9 sigma apart
    files = sorted(glob.glob(str(good / "data_collection" / "window*.ptRes")))
    ours = neighbour_overlaps(files)
    browser = [o for _a, _b, o in USRun(str(good)).overlap_1d()[0]]
    assert len(ours) == 19 and max(abs(a - b) for a, b in zip(ours, browser)) < 1e-12
    assert min(ours) > 0.2
    assert max(neighbour_overlaps(sorted(glob.glob(str(stiff / "data_collection" / "window*.ptRes"))))) < 0.03


_OPT_LOG = """
----------------------------------------------------------------------------------------------------------------
  Iteration       Function          RMS Gradient        Max. |Grad.|          RMS Disp.         Max. |Disp.|
----------------------------------------------------------------------------------------------------------------
     0     I         52.31646236        278.55447516       1485.02299486          0.01230915          0.06562224
    50   L1s       -108.07687135          0.89258957          2.44377627          0.00094669          0.00204934
   100   L1s       -108.17965437          0.09488795          1.10177029          0.00054205          0.00172325
----------------------------------------------------------------------------------------------------------------

Minimization converged.

------------------------------- Summary of Energy Terms ------------------------------
Potential Energy        =        -108.1797  RMS Gradient            =             None
Harmonic Angle          =           3.2471  Harmonic Bond           =           0.9751
--------------------------------------------------------------------------------------
"""

_MD_LOG = """
        Time            Total Energy       Kinetic Energy     Potential Energy       Temperature
----------------------------------------------------------------------------------------------------
          0.00000000       -256.97934707         74.83025420       -331.80960127        300.00000000
          0.10000000       -244.65746205         49.26691793       -293.92437999        200.00000000
          0.20000000       -214.45923518         77.14715594       -291.60639112        310.00000000
          0.30000000       -209.00959673         63.59799749       -272.60759422        250.00000000
----------------------------------------------------------------------------------------------------
"""


def test_log_analysis_kinds(tmp_path):
    from util.log_analysis import analyze, analyze_xtb
    opt = tmp_path / "opt.log"
    opt.write_text(_OPT_LOG)
    r = analyze(str(opt))
    s = dict(r["summary"])
    assert r["kind"] == "optimization" and s["status"] == "Minimization converged." and s["iterations"] == "100"
    assert s["final energy"].startswith("-108.1797") and len(r["panels"]) == 3     # . energy, gradient, terms
    md = tmp_path / "md" / "output.log"
    md.parent.mkdir()
    md.write_text(_MD_LOG)
    r = analyze(str(md), "Molecular_Dynamics")
    assert r["kind"] == "md" and [p["title"] for p in r["panels"]] == ["Energies", "Temperature"]
    assert dict(r["summary"])["temperature"].startswith("265.00")                     # . 4 points: none dropped
    for i, (e, g) in enumerate(((-10.0, 5.0), (-10.001, 5.5))):
        (md.parent / "frame{}.xtb.out".format(i)).write_text(
            ":: HOMO-LUMO gap  {} eV ::\n| TOTAL ENERGY  {} Eh |\n".format(g, e))
    assert analyze(str(md), "Molecular_Dynamics")["qc"]["program"] == "xTB"
    x = analyze_xtb(sorted(str(p) for p in md.parent.glob("*.xtb.out")))
    assert dict(x["summary"])["HOMO–LUMO gap"] == "5.00 – 5.50 eV"
    bad = tmp_path / "unknown.log"
    bad.write_text("nothing here\n")
    assert analyze(str(bad))["kind"] is None
