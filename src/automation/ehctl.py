#!/usr/bin/env python3
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  ehctl -- command-line client of the EasyHybrid control server.
#
#    python3 src/automation/ehctl.py launch                # start EasyHybrid --control
#    python3 src/automation/ehctl.py state
#    python3 src/automation/ehctl.py help                  # commands
#    python3 src/automation/ehctl.py help scan             # one command
#    python3 src/automation/ehctl.py load path=examples/tutorials/01_bala/bALA_not_opt.mol
#    python3 src/automation/ehctl.py atoms obj=0 resn=ALA limit=5
#    python3 src/automation/ehctl.py optimize obj=0 rms_gradient=0.05
#    python3 src/automation/ehctl.py wait sys=0
#    python3 src/automation/ehctl.py scan --json '{"rc1": {"atoms": [4, 16], "step": 0.1, "nsteps": 10}}'
#    python3 src/automation/ehctl.py terminal command="show sticks obj=0"
#    python3 src/automation/ehctl.py screenshot path=/tmp/view.png
#
#  key=value arguments are read as JSON when possible (numbers, true/false,
#  null, [lists], {dicts}); otherwise as text. Output: JSON.

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from automation.client import EasyHybridClient, EasyHybridError, launch   # noqa: E402


USAGE = """usage: ehctl <command> [key=value ...] [--json '{...}'] [--timeout s]

  ehctl launch                 start EasyHybrid with the control server
  ehctl help [command]         list the commands / arguments of one
  ehctl state                  systems, objects, jobs, dialogs
  ehctl wait sys=0 [job=3]     wait for a job

Values are read as JSON when possible (numbers, true/false, [lists], {dicts})."""


def _value(text):
    try:
        return json.loads(text)
    except ValueError:
        return text


def _help(client, name=None):
    schemas = client.call("describe")
    schemas += [{"name": "wait", "summary": "Waits for a job to finish (sys=, job=, timeout= s).",
                 "description": "", "inputSchema": {"properties": {}}},
                {"name": "result", "summary": "Result of a pending request (request=r7).",
                 "description": "", "inputSchema": {"properties": {}}}]
    if name is None:
        width = max(len(s["name"]) for s in schemas)
        lines = ["Commands (ehctl help <command> for its arguments):"]
        lines += ["  {:<{w}}  {}".format(s["name"], s["summary"], w=width) for s in sorted(schemas, key=lambda s: s["name"])]
        return "\n".join(lines)
    for s in schemas:
        if s["name"] == name:
            lines = [name, "-" * len(name), s["description"] or s["summary"], "", "Arguments:"]
            props = s["inputSchema"].get("properties", {})
            required = set(s["inputSchema"].get("required", []))
            for arg, p in props.items():
                default = "" if arg in required else " = {}".format(json.dumps(p.get("default")))
                lines.append("  {}{}{}  {}".format(arg, default, " (required)" if arg in required else "",
                                                   p.get("description", "")))
            return "\n".join(lines)
    return "No command '{}'.".format(name)


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(USAGE)
        return 0
    command, rest = argv[0], argv[1:]
    if command == "launch":
        process = launch(rest)
        print(json.dumps({"started": process is not None, "pid": process.pid if process else None}))
        return 0
    params, timeout = {}, 60
    i = 0
    while i < len(rest):
        arg = rest[i]
        if arg == "--json":
            params.update(json.loads(rest[i + 1]))
            i += 2
            continue
        if arg == "--timeout":
            timeout = float(rest[i + 1])
            i += 2
            continue
        if "=" not in arg:
            params.setdefault("_args", []).append(arg)
            i += 1
            continue
        key, value = arg.split("=", 1)
        params[key] = _value(value)
        i += 1
    client = EasyHybridClient()
    try:
        if command == "help":
            print(_help(client, (params.get("_args") or [None])[0]))
            return 0
        params.pop("_args", None)
        if command == "wait":
            params.setdefault("timeout", 3600)
            timeout = float(params["timeout"]) + 30
        reply = client.request(command, params, timeout)
    except EasyHybridError as error:
        print(json.dumps({"ok": False, "error": str(error)}), file=sys.stderr)
        return 2
    finally:
        client.close()
    print(json.dumps(reply, indent=2, default=str))
    return 0 if reply.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
