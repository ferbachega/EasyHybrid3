#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Client of the control server (standard library only; usable from any
#  Python, without GTK or pDynamo):
#
#      from automation.client import EasyHybridClient
#      eh = EasyHybridClient()            # ~/.easyhybrid/control.sock
#      eh.call("state")
#      job = eh.call("optimize", obj=0)
#      eh.wait(sys=job["system"], job=job["job"])
#
#  launch() starts EasyHybrid with the control server and waits for it.

import itertools
import json
import os
import socket
import subprocess
import sys
import time

CONTROL_DIR = os.path.join(os.path.expanduser("~"), ".easyhybrid")
DEFAULT_SOCKET = os.path.join(CONTROL_DIR, "control.sock")
EASYHYBRID_HOME = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class EasyHybridError(Exception):
    """ The server replied ok=false (message = its error). """

    def __init__(self, message, reply=None):
        super().__init__(message)
        self.reply = reply or {}


def socket_path():
    """ The first existing of: EASYHYBRID_SOCKET, the socket of the last
        EasyHybrid started with the control server (~/.easyhybrid/control.json;
        the server moves a too-long path to a short one), the default. """
    candidates = [os.environ.get("EASYHYBRID_SOCKET")]
    try:
        with open(os.path.join(CONTROL_DIR, "control.json")) as f:
            candidates.append(json.load(f).get("socket"))
    except (OSError, ValueError):
        pass
    candidates.append(DEFAULT_SOCKET)
    for path in candidates:
        if path and os.path.exists(path):
            return path
    return candidates[0] or DEFAULT_SOCKET


class EasyHybridClient:
    def __init__(self, path=None, timeout=None):
        self.path = path or socket_path()
        self.timeout = timeout
        self._sock = None
        self._reader = None
        self._ids = itertools.count(1)

    def connect(self):
        if self._sock is None:
            sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            sock.settimeout(self.timeout)
            try:
                sock.connect(self.path)
            except OSError as error:
                sock.close()
                raise EasyHybridError("EasyHybrid is not running with the control server ({}: {}). "
                                      "Start it with: python3 easyhybrid.py --control".format(self.path, error))
            self._sock = sock
            self._reader = sock.makefile("r", encoding="utf-8")
        return self

    def close(self):
        if self._sock is not None:
            try:
                self._reader.close()
                self._sock.close()
            finally:
                self._sock = self._reader = None

    def __enter__(self):
        return self.connect()

    def __exit__(self, *exc):
        self.close()

    def request(self, method, params=None, timeout=60):
        """ Sends one request; returns the full reply dict. """
        self.connect()
        message = {"id": next(self._ids), "method": method, "params": params or {}, "timeout": timeout}
        self._sock.sendall((json.dumps(message) + "\n").encode("utf-8"))
        line = self._reader.readline()
        if not line:
            self.close()
            raise EasyHybridError("Connection closed by EasyHybrid.")
        return json.loads(line)

    def call(self, _command, _timeout=60, **params):
        """ Runs a command; returns its result (raises EasyHybridError on
            failure). A pending reply (dialog open) is returned as is.
            _timeout: seconds the server waits before replying "pending". """
        reply = self.request(_command, params, _timeout)
        if not reply.get("ok"):
            raise EasyHybridError(reply.get("error", "error"), reply)
        if reply.get("pending"):
            return reply
        return reply.get("result")

    def wait(self, sys=None, job=None, timeout=3600):
        """ Waits for a job to finish; returns its status. """
        return self.call("wait", _timeout=timeout + 30, sys=sys, job=job, timeout=timeout)


def launch(args=(), wait=90, python=None, log=None):
    """ Starts EasyHybrid with the control server; returns the Popen once
        the socket answers (raises EasyHybridError after `wait` seconds). """
    path = socket_path()
    probe = EasyHybridClient(path, timeout=2)
    try:
        probe.call("ping")
        probe.close()
        return None                     # already running
    except (EasyHybridError, OSError):
        probe.close()
    env = dict(os.environ)
    env["EASYHYBRID_SOCKET"] = env.get("EASYHYBRID_SOCKET") or DEFAULT_SOCKET
    out = open(log, "a") if log else subprocess.DEVNULL
    process = subprocess.Popen([python or sys.executable, os.path.join(EASYHYBRID_HOME, "easyhybrid.py"),
                                "--control"] + list(args),
                               cwd=EASYHYBRID_HOME, env=env, stdout=out, stderr=subprocess.STDOUT,
                               start_new_session=True)
    deadline = time.time() + wait
    while time.time() < deadline:
        if process.poll() is not None:
            raise EasyHybridError("EasyHybrid exited during start (code {}).".format(process.returncode))
        try:
            client = EasyHybridClient(timeout=2)
            client.call("ping")
            client.close()
            return process
        except (EasyHybridError, OSError):
            time.sleep(1.0)
    raise EasyHybridError("EasyHybrid did not open the control socket in {} s.".format(wait))
