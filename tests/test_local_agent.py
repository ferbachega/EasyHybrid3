"""
Testes do assistente local (src/automation/local_agent.py) com um servidor
HTTP falso que imita o Ollama e o formato OpenAI: laco de ferramentas,
confirmacao de acoes, ferramentas proibidas e a lista de modelos.
"""
import json
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

import pytest

from automation import local_agent


class _Server:
    def __init__(self, replies, models=None):
        self.replies, self.requests, self.models = list(replies), [], models or {}
        server = self

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def _send(self, payload):
                body = json.dumps(payload).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def do_GET(self):
                self._send(server.models.get(self.path, {}))

            def do_POST(self):
                data = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
                server.requests.append((self.path, data))
                self._send(server.replies.pop(0))

        self.httpd = HTTPServer(("127.0.0.1", 0), Handler)
        self.url = "http://127.0.0.1:{}".format(self.httpd.server_port)
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()

    def close(self):
        self.httpd.shutdown()


def _ollama_call(name, args):
    return {"message": {"role": "assistant", "content": "", "tool_calls": [{"function": {"name": name, "arguments": args}}]}}


def _ollama_answer(text):
    return {"message": {"role": "assistant", "content": text}}


def test_tool_loop_ollama():
    server = _Server([_ollama_call("pes", {"logfile": "/x/output.log"}), _ollama_answer("barrier 8.3 kJ/mol")])
    calls, events = [], []
    agent = local_agent.LocalAgent("ollama", server.url, "m", execute=lambda n, a: calls.append((n, a)) or {"barrier": 8.3},
                                   on_event=lambda k, d: events.append(k))
    try:
        assert agent.ask("barrier?") == "barrier 8.3 kJ/mol"
    finally:
        server.close()
    assert calls == [("pes", {"logfile": "/x/output.log"})]
    path, first = server.requests[0]
    assert path == "/api/chat" and first["stream"] is False and first["options"]["num_ctx"] == 16384
    names = {t["function"]["name"] for t in first["tools"]}
    assert {"pes", "state", "scan", "wait", "help"} <= names and not names & {"python", "simulate", "screenshot"}
    tool_msg = server.requests[1][1]["messages"][-1]
    assert tool_msg["role"] == "tool" and tool_msg["tool_name"] == "pes" and json.loads(tool_msg["content"])["ok"]
    assert events.count("tool_call") == 1 and events[-1] == "answer"


def test_actions_need_confirmation_and_excluded_tools():
    server = _Server([_ollama_call("optimize", {"obj": 0}), _ollama_call("python", {"code": "1"}), _ollama_answer("ok")])
    executed = []
    agent = local_agent.LocalAgent("ollama", server.url, "m", execute=lambda n, a: executed.append(n),
                                   confirm=lambda n, a: False)
    try:
        agent.ask("optimize it")
    finally:
        server.close()
    assert executed == []
    results = [json.loads(m["content"]) for _p, r in server.requests[1:] for m in r["messages"][-1:]]
    assert "did not allow" in results[0]["error"] and "not available" in results[1]["error"]


def test_openai_format():
    replies = [{"choices": [{"message": {"role": "assistant", "content": None, "tool_calls": [
                   {"id": "c1", "type": "function", "function": {"name": "state", "arguments": "{}"}}]}}]},
               {"choices": [{"message": {"role": "assistant", "content": "2 systems"}}]}]
    server = _Server(replies, models={"/v1/models": {"data": [{"id": "qwen2.5-7b"}]}})
    agent = local_agent.LocalAgent("openai", server.url, "qwen", execute=lambda n, a: {"systems": 2})
    try:
        assert agent.ask("how many systems?") == "2 systems"
        assert local_agent.list_models("openai", server.url) == [("qwen2.5-7b", "")]
    finally:
        server.close()
    assert server.requests[0][0] == "/v1/chat/completions"
    tool_msg = server.requests[1][1]["messages"][-1]
    assert tool_msg == {"role": "tool", "tool_call_id": "c1", "content": tool_msg["content"]}


def test_ollama_model_list_and_unreachable():
    server = _Server([], models={"/api/tags": {"models": [
        {"name": "gpt-oss:20b", "details": {"parameter_size": "20.9B"}, "capabilities": ["completion", "tools"]},
        {"name": "tiny", "details": {"parameter_size": "1B"}, "capabilities": ["completion"]}]}})
    try:
        assert local_agent.list_models("ollama", server.url) == [("gpt-oss:20b", "20.9B, tools"), ("tiny", "1B, no tool support")]
    finally:
        server.close()
    with pytest.raises(local_agent.AgentError, match="Cannot reach"):
        local_agent.list_models("ollama", "http://127.0.0.1:9")
