#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Local AI assistant (2026-10-28): a chat with a model served on this
#  computer (Ollama, or any server with the OpenAI chat API: LM Studio,
#  llama.cpp, vLLM) that operates EasyHybrid through the same commands as the
#  MCP server (automation/api.py). Standard library only (urllib).
#
#  The agent runs in a worker thread. It asks `execute(name, args)` to run a
#  command (the window runs it in the GTK thread when needed) and
#  `confirm(name, args)` before an action (start a simulation, load a file,
#  change the QC region...). Commands that read only run without asking.
#  Some commands are never offered to a local model (python, simulate, remote
#  machine setup, screenshots: a text model cannot see images).

import json
import time
import urllib.error
import urllib.request

from automation import settings as ai_settings

# . what a local model may do
READ_ONLY = {"ping", "state", "atoms", "measure", "energy", "pes", "job_summary", "us_check", "jobs", "job_status",
             "events", "dialogs", "windows", "remotes", "remote_jobs", "remote_partitions", "set_visible", "set_frame",
             "us_browser"}
ACTIONS = {"load", "new_system", "save_session", "activate", "set_qc", "fix_atoms", "optimize", "scan", "md", "neb",
           "umbrella", "wham", "import_result", "terminal", "run_on", "remote_fetch", "remote_cancel", "dialog_respond",
           "result"}
EXCLUDED = {"python", "simulate", "remote_save", "remote_test", "remote_detect", "screenshot", "launch"}

PROVIDERS = (("ollama", "Ollama"), ("openai", "OpenAI-compatible server (LM Studio, llama.cpp, vLLM)"))
DEFAULT_URLS = {"ollama": "http://localhost:11434", "openai": "http://localhost:1234"}

SYSTEM = """ALWAYS answer in the same language as the user's last message (English if they write English,
Portuguese if they write Portuguese, and so on). Never switch to another language.
You are the assistant of EasyHybrid, a program for QC/MM and molecular simulations (pDynamo3).
You act through the tools, which operate the EasyHybrid window the user is looking at.
- When the user gives the path of a log file, read it directly (job_summary or pes with logfile=);
  do not look up the session first.
- Start with `state` when you need ids: systems (sys=) and objects (obj=). Atom indices are 0-based.
- Simulations (optimize, scan, md, neb, umbrella) start a job; with wait=true the call returns its
  result. Otherwise call `wait`, then `import_result` and `pes` / `job_summary`.
- Umbrella sampling: check `us_check` before trusting a PMF ("fail" means it is not a result).
- Energies are in kJ/mol, distances in Angstrom. Write "QC/MM", never "QM/MM".
- Before an action the user is asked to confirm it; if it is refused, do not repeat it, ask what to do.
- Be brief. Give the numbers. Answer in the language of the user.
- You cannot see images; ask the user to look at the 3D view when needed."""


class AgentError(Exception):
    pass


def _post(url, payload, timeout):
    req = urllib.request.Request(url, data=json.dumps(payload).encode("utf-8"),
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as error:
        raise AgentError("{} answered {}: {}".format(url, error.code, error.read().decode("utf-8", "replace")[:300]))
    except (urllib.error.URLError, OSError) as error:
        raise AgentError("Cannot reach {} ({}). Is the server running?".format(url, error))


def _get(url, timeout=10):
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, OSError, ValueError) as error:
        raise AgentError("Cannot reach {} ({}).".format(url, error))


def list_models(provider, url):
    """ [(name, note)] of the models the server offers (Ollama: with "tools" when it supports them). """
    url = url.rstrip("/")
    if provider == "ollama":
        data = _get(url + "/api/tags")
        out = []
        for m in data.get("models", []):
            caps = m.get("capabilities") or (m.get("details") or {}).get("capabilities") or []
            size = (m.get("details") or {}).get("parameter_size") or ""
            out.append((m["name"], "{}{}".format(size, ", tools" if "tools" in caps else (", no tool support" if caps else ""))))
        return out
    data = _get(url + "/v1/models")
    return [(m.get("id"), "") for m in data.get("data", [])]


def tool_specs(names=None):
    """ The commands a local model gets, in the OpenAI/Ollama "tools" format
        (short descriptions: local models have small contexts). """
    from automation.api import command_schemas
    from automation.mcp_server import _short_schema
    specs = []
    for s in command_schemas():
        name = s["name"]
        if name in EXCLUDED or (names is not None and name not in names):
            continue
        if name not in READ_ONLY and name not in ACTIONS:
            continue
        schema = _short_schema(s["inputSchema"])
        # . local models get summaries only (small contexts): no full data tables, no image sizes
        for arg in ("detail", "max_width"):
            schema.get("properties", {}).pop(arg, None)
        specs.append({"type": "function", "function": {"name": name, "description": s["summary"],
                                                        "parameters": schema}})
    specs.append({"type": "function", "function": {
        "name": "wait", "description": "Waits until a job finishes; returns its status.",
        "parameters": {"type": "object", "properties": {"sys": {"type": "integer"}, "job": {"type": "integer"},
                                                        "timeout": {"type": "number", "default": 600}}}}})
    specs.append({"type": "function", "function": {
        "name": "help", "description": "Full description and arguments of a tool.",
        "parameters": {"type": "object", "properties": {"tool": {"type": "string"}}, "required": ["tool"]}}})
    return specs


class LocalAgent:
    """ One conversation. ask() runs the tool loop and returns the final text. """

    def __init__(self, provider="ollama", url=None, model=None, execute=None, confirm=None, on_event=None,
                 num_ctx=16384, max_steps=12, timeout=600, settings=None):
        self.provider = provider
        self.url = (url or DEFAULT_URLS.get(provider, "")).rstrip("/")
        self.model = model
        self.execute = execute
        self.confirm = confirm
        self.on_event = on_event or (lambda kind, data: None)
        self.num_ctx = int(num_ctx)
        self.max_steps = int(max_steps)
        self.timeout = timeout
        self.settings = settings if settings is not None else ai_settings.load()
        self.tools = tool_specs()
        self.messages = [{"role": "system", "content": SYSTEM}]
        self.stopped = False

    # ------------------------------------------------------------------ server
    def _chat(self):
        if self.provider == "ollama":
            data = _post(self.url + "/api/chat", {"model": self.model, "messages": self.messages, "tools": self.tools,
                                                  "stream": False, "options": {"num_ctx": self.num_ctx}}, self.timeout)
            message = data.get("message") or {}
            calls = [{"id": None, "name": c["function"]["name"], "arguments": c["function"].get("arguments") or {}}
                     for c in message.get("tool_calls") or []]
            self.messages.append({k: v for k, v in message.items() if k in ("role", "content", "tool_calls", "thinking")})
            return message.get("content") or "", calls
        data = _post(self.url + "/v1/chat/completions", {"model": self.model, "messages": self.messages,
                                                         "tools": self.tools}, self.timeout)
        message = (data.get("choices") or [{}])[0].get("message") or {}
        calls = []
        for c in message.get("tool_calls") or []:
            args = c["function"].get("arguments") or "{}"
            try:
                args = json.loads(args) if isinstance(args, str) else args
            except ValueError:
                args = {}
            calls.append({"id": c.get("id"), "name": c["function"]["name"], "arguments": args})
        self.messages.append({"role": "assistant", "content": message.get("content") or "",
                              **({"tool_calls": message["tool_calls"]} if message.get("tool_calls") else {})})
        return message.get("content") or "", calls

    def _tool_message(self, call, content):
        if self.provider == "ollama":
            return {"role": "tool", "content": content, "tool_name": call["name"]}
        return {"role": "tool", "tool_call_id": call["id"], "content": content}

    # ------------------------------------------------------------------ tools
    def _run_tool(self, call):
        name, args = call["name"], dict(call["arguments"] or {})
        args.pop("detail", None)                     # . summaries only for a local model
        if name in EXCLUDED or (name not in READ_ONLY and name not in ACTIONS and name not in ("wait", "help")):
            return {"ok": False, "error": "The tool '{}' is not available to the local assistant.".format(name)}
        if name == "help":
            from automation.api import command_schemas
            for s in command_schemas():
                if s["name"] == args.get("tool"):
                    return {"ok": True, "result": {"description": s["description"], "arguments": s["inputSchema"]}}
            return {"ok": False, "error": "No tool named {}".format(args.get("tool"))}
        if name in ACTIONS and self.confirm is not None and not self.confirm(name, args):
            return {"ok": False, "error": "The user did not allow this action."}
        if name == "wait":
            return self._wait(args)
        try:
            return {"ok": True, "result": self.execute(name, args)}
        except Exception as error:
            return {"ok": False, "error": "{}: {}".format(type(error).__name__, error)}

    def _wait(self, args):
        deadline = time.time() + float(args.get("timeout") or 600)
        while not self.stopped:
            try:
                info = self.execute("job_status", {k: args[k] for k in ("sys", "job") if k in args})
            except Exception as error:
                return {"ok": False, "error": str(error)}
            if info.get("status") not in (None, "Queued", "Running", "Running..."):
                return {"ok": True, "result": info}
            if time.time() > deadline:
                info["timeout"] = True
                return {"ok": True, "result": info}
            time.sleep(2.0)
        return {"ok": False, "error": "stopped"}

    # ------------------------------------------------------------------ loop
    def ask(self, text):
        """ Adds the user's message and runs model -> tools -> model until the
            model answers without tools (or max_steps). Returns the answer. """
        if not self.model:
            raise AgentError("Choose a model first.")
        self.stopped = False
        self.messages.append({"role": "user", "content": text})
        for _step in range(self.max_steps):
            if self.stopped:
                return "(stopped)"
            started = time.time()
            content, calls = self._chat()
            self.on_event("model", {"seconds": time.time() - started, "tools": [c["name"] for c in calls]})
            if not calls:
                self.on_event("answer", {"text": content})
                return content
            if content.strip():
                self.on_event("thought", {"text": content})
            for call in calls:
                self.on_event("tool_call", {"name": call["name"], "arguments": call["arguments"]})
                reply = self._run_tool(call)
                text_reply = json.dumps(reply, separators=(",", ":"), default=str)
                if len(text_reply) > 8000:
                    text_reply = text_reply[:8000] + "...(truncated)"
                self.on_event("tool_result", {"name": call["name"], "ok": reply.get("ok"),
                                              "text": text_reply[:600]})
                self.messages.append(self._tool_message(call, text_reply))
        note = "(stopped after {} steps without a final answer)".format(self.max_steps)
        self.on_event("answer", {"text": note})
        return note
