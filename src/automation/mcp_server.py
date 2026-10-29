#!/usr/bin/env python3
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Model Context Protocol server (stdio) for AI assistants (2026-10-28).
#  Standard library only. Each tool forwards to the control server of a
#  running EasyHybrid (the commands of automation/api.py), plus:
#    launch      start EasyHybrid with the control server
#    wait        wait for a job
#    result      result of a pending request (a dialog was open)
#  screenshot also returns the image, so the assistant sees it.
#
#  Claude Code:   claude mcp add easyhybrid -- python3 /path/to/EasyHybrid/src/automation/mcp_server.py
#  Other clients: command "python3", args ["/path/to/EasyHybrid/src/automation/mcp_server.py"]

import base64
import json
import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automation import settings as ai_settings                                # noqa: E402
from automation.api import command_schemas                                    # noqa: E402
from automation.client import EasyHybridClient, EasyHybridError, launch       # noqa: E402

PROTOCOL_VERSION = "2025-06-18"

INSTRUCTIONS = """EasyHybrid is a GUI for QC/MM and molecular simulations (pDynamo3).
These tools drive a running EasyHybrid window; call `launch` first if it is not running.
Start with `state`: it lists systems (sys=id), objects (obj=id: coordinate sets and trajectories),
running jobs and open dialogs. Atom indices are 0-based pDynamo indices (see `atoms`).
Simulations (optimize, scan, md, neb, simulate) return a job at once: call `wait`, then
`import_result` to get the trajectory as a new object, and `pes` for the energies as numbers.
Umbrella sampling: `umbrella` (windows from a finished scan job), `wait`, then `us_check` (do not
trust the PMF when its verdict is "fail") and `wham` for the free energy profile.
Messages that dialogs would show come back in `events`; a question dialog makes a reply
"pending": read it with `dialogs`, answer with `dialog_respond`, then call `result`.
Use `screenshot` to look at the 3D view or a window, and `terminal` for display commands
(show/hide/color/center/select; `terminal command="help"`). Energies are in kJ/mol, distances in A."""

EXTRA_TOOLS = [
    {"name": "launch", "description": "Starts EasyHybrid (GUI) with the control server, if it is not running.",
     "inputSchema": {"type": "object", "properties": {
         "file": {"type": "string", "description": "optional .easy/.pkl/.pdb to open at start"}}}},
    {"name": "wait", "description": "Waits until a job finishes (or timeout s); returns its status.",
     "inputSchema": {"type": "object", "properties": {
         "sys": {"type": ["integer", "string"], "description": "system id (default: the active system)"},
         "job": {"type": "integer", "description": "job number (default: the last job)"},
         "timeout": {"type": "number", "default": 1800}}}},
    {"name": "help", "description": "Full description and arguments of an EasyHybrid tool.",
     "inputSchema": {"type": "object", "properties": {"tool": {"type": "string"}}, "required": ["tool"]}},
    {"name": "result", "description": "Result of a pending request (reply had \"pending\": \"r7\").",
     "inputSchema": {"type": "object", "properties": {"request": {"type": "string"},
                                                     "timeout": {"type": "number", "default": 30}},
                     "required": ["request"]}},
]


def _short_schema(schema, limit=70):
    """ short_tools: descriptions only for the required arguments (cut to one
        line); the conventions of sys/obj/frame/atoms are in the instructions
        and the rest comes from the help tool. """
    import copy
    schema = copy.deepcopy(schema)
    required = set(schema.get("required", []))
    for name, prop in schema.get("properties", {}).items():
        text = prop.pop("description", None)
        if text and name in required:
            prop["description"] = text if len(text) <= limit else text[:limit - 3].rstrip() + "..."
    return schema


def tool_list(settings=None):
    settings = settings or ai_settings.load()
    short = settings.get("short_tools")
    tools = []
    for s in command_schemas():
        if short:
            tools.append({"name": s["name"], "description": s["summary"], "inputSchema": _short_schema(s["inputSchema"])})
        else:
            tools.append({"name": s["name"], "description": s["description"] or s["summary"],
                          "inputSchema": s["inputSchema"]})
    return tools + EXTRA_TOOLS


def instructions(settings=None):
    """ INSTRUCTIONS plus what the user chose in Edit > AI Control Server. """
    settings = settings or ai_settings.load()
    text = INSTRUCTIONS
    if settings.get("summaries"):
        text += ("\nThe user chose economical token use: replies are summaries; ask detail=\"full\" only when "
                 "the data arrays are really needed. Prefer numbers to screenshots.")
    if settings.get("auto_wait"):
        text += ("\nSimulations wait for their job and return the imported result and its summary in the same "
                 "call (wait=false to start without waiting).")
    if settings.get("short_tools"):
        text += "\nTool descriptions are short: call `help` with a tool name for its full description."
    model = settings.get("operator_model", "session")
    if model != "session":
        text += ("\nThe user chose Claude {} to operate EasyHybrid: delegate the execution steps (loading, running "
                 "and waiting for simulations, collecting numbers) to the `easyhybrid-operator` subagent, and keep "
                 "planning and interpretation in this session.".format(model.capitalize()))
    return text


class MCPServer:
    def __init__(self):
        self.client = None

    def _client(self):
        if self.client is None:
            self.client = EasyHybridClient().connect()
        return self.client

    def _forward(self, method, params, timeout=120):
        try:
            return self._client().request(method, params, timeout)
        except (EasyHybridError, OSError):
            # . reconnect once (EasyHybrid restarted)
            if self.client is not None:
                self.client.close()
            self.client = None
            return self._client().request(method, params, timeout)

    def call_tool(self, name, args):
        args = dict(args or {})
        if name == "help":
            for s in command_schemas():
                if s["name"] == args.get("tool"):
                    return [{"type": "text", "text": json.dumps({"name": s["name"], "description": s["description"],
                                                                 "arguments": s["inputSchema"]}, indent=1)}], False
            return [{"type": "text", "text": "No tool named {}.".format(args.get("tool"))}], True
        if name == "launch":
            extra = [args["file"]] if args.get("file") else []
            process = launch(extra)
            return [{"type": "text", "text": json.dumps({"started": process is not None,
                                                          "note": "already running" if process is None else "ready"})}], False
        timeout = 120
        if name == "wait":
            args.setdefault("timeout", 1800)
            timeout = float(args["timeout"]) + 30
        if name in ("optimize", "scan", "md", "neb", "umbrella"):
            timeout = 3600                    # . auto-wait may keep the call open
        reply = self._forward(name, args, timeout)
        if reply.get("format") == "compact":
            reply.pop("format", None)
            text = json.dumps(reply, separators=(",", ":"), default=str)
        else:
            text = json.dumps(reply, indent=1, default=str)
        content = [{"type": "text", "text": text}]
        if name == "screenshot" and reply.get("ok"):
            path = (reply.get("result") or {}).get("path")
            if path and os.path.isfile(path) and os.path.getsize(path) < 8 * 1024 * 1024:
                with open(path, "rb") as f:
                    content.append({"type": "image", "mimeType": "image/png",
                                    "data": base64.b64encode(f.read()).decode("ascii")})
        return content, not reply.get("ok", False)

    def handle(self, message):
        method = message.get("method")
        mid = message.get("id")
        params = message.get("params") or {}
        if mid is None:
            return None                                   # notification
        if method == "initialize":
            result = {"protocolVersion": params.get("protocolVersion") or PROTOCOL_VERSION,
                      "capabilities": {"tools": {"listChanged": False}},
                      "serverInfo": {"name": "easyhybrid", "version": "3.2"},
                      "instructions": instructions()}
        elif method == "ping":
            result = {}
        elif method == "tools/list":
            result = {"tools": tool_list()}
        elif method == "tools/call":
            try:
                content, is_error = self.call_tool(params.get("name"), params.get("arguments"))
            except EasyHybridError as error:
                content, is_error = [{"type": "text", "text": str(error)}], True
            except Exception as error:
                content = [{"type": "text", "text": "{}: {}\n{}".format(type(error).__name__, error,
                                                                       traceback.format_exc())}]
                is_error = True
            result = {"content": content, "isError": is_error}
        else:
            return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32601, "message": "Method not found: {}".format(method)}}
        return {"jsonrpc": "2.0", "id": mid, "result": result}

    def serve(self, stdin=sys.stdin, stdout=sys.stdout):
        for line in stdin:
            line = line.strip()
            if not line:
                continue
            try:
                message = json.loads(line)
            except ValueError:
                reply = {"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": "Parse error"}}
            else:
                reply = self.handle(message)
            if reply is not None:
                stdout.write(json.dumps(reply) + "\n")
                stdout.flush()


if __name__ == "__main__":
    MCPServer().serve()
