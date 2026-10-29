#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Control server (2026-10-28): a Unix socket inside the running EasyHybrid.
#
#  Protocol: one JSON object per line, both ways.
#    request  {"id": 1, "method": "state", "params": {...}, "timeout": 60}
#    reply    {"id": 1, "ok": true, "result": ..., "events": [...]}
#             {"id": 1, "ok": false, "error": "...", "traceback": "...", "events": [...]}
#             {"id": 1, "ok": true, "pending": "r7", "dialogs": [...], "events": [...]}
#
#  Commands are the methods of AutomationAPI and run in the GTK main thread
#  (GLib.idle_add). A command that opens a modal dialog does not block the
#  client: after a moment the reply is "pending" with the dialog; answer it
#  with dialog_respond and collect the result with {"method": "result",
#  "params": {"request": "r7"}}. Handled here, in the client thread:
#    describe   the command list with argument schemas
#    wait       waits for a job: {"sys": 0, "job": 3, "timeout": 600}
#    result     result of a pending request
#
#  The socket is created with mode 0600 in ~/.easyhybrid (0700): only the
#  same user can connect. Started by `easyhybrid.py --control`, by the
#  EASYHYBRID_CONTROL=1 environment variable or by the terminal command
#  `control on`.

import itertools
import json
import os
import socket
import tempfile
import threading
import time
import traceback

from automation import events
from automation import settings as ai_settings
from automation.api import APIError, AutomationAPI, _jsonable, command_schemas

# . commands whose job the server can wait for and summarize in the same call
SIMULATIONS = ("optimize", "scan", "md", "neb", "umbrella")

CONTROL_DIR = os.path.join(os.path.expanduser("~"), ".easyhybrid")
DEFAULT_SOCKET = os.path.join(CONTROL_DIR, "control.sock")
INFO_FILE = os.path.join(CONTROL_DIR, "control.json")


def default_socket_path():
    return os.environ.get("EASYHYBRID_SOCKET") or DEFAULT_SOCKET


def _strip_empty(value):
    """ Compact replies: None values and empty lists/dicts left out. """
    if isinstance(value, dict):
        return {k: _strip_empty(v) for k, v in value.items() if v is not None and v != [] and v != {}}
    if isinstance(value, list):
        return [_strip_empty(v) for v in value]
    return value


def _alive(path):
    s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    s.settimeout(1.0)
    try:
        s.connect(path)
        return True
    except OSError:
        return False
    finally:
        s.close()


class ControlServer:
    """ Unix-socket server that runs AutomationAPI commands in the GUI. """

    def __init__(self, main, path=None, dialogs_off=None, exit_handlers=True, settings=None):
        self.main = main
        self.exit_handlers = exit_handlers
        self.settings = dict(settings) if settings is not None else ai_settings.load()
        self.api = AutomationAPI(main, settings=self.settings)
        self.path = path or default_socket_path()
        # . None: from the AI settings ("dialogs")
        self.dialogs_off = (not self.settings.get("dialogs")) if dialogs_off is None else dialogs_off
        self._sock = None
        self._thread = None
        self._running = False
        self._pending = {}
        self._rid = itertools.count(1)

    # ---------------------------------------------------------------- life
    def start(self):
        if len(os.fsencode(self.path)) > 100:
            # . AF_UNIX paths are limited to ~107 bytes: use a short one (the
            #   clients find it in ~/.easyhybrid/control.json)
            self.path = os.path.join(tempfile.gettempdir(), "easyhybrid-{}".format(os.getuid()),
                                     "control-{}.sock".format(os.getpid()))
        os.makedirs(os.path.dirname(self.path), mode=0o700, exist_ok=True)
        os.makedirs(CONTROL_DIR, mode=0o700, exist_ok=True)
        for folder in {os.path.dirname(self.path), CONTROL_DIR,
                       os.path.join(tempfile.gettempdir(), "easyhybrid-{}".format(os.getuid()))}:
            if os.path.isdir(folder):
                self._remove_stale_sockets(folder)
        if os.path.exists(self.path):
            if _alive(self.path):
                # . another EasyHybrid already listens there: use a per-process socket
                base, ext = os.path.splitext(self.path)
                self.path = "{}-{}{}".format(base, os.getpid(), ext or ".sock")
            else:
                os.unlink(self.path)
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        old_umask = os.umask(0o177)
        try:
            sock.bind(self.path)
        finally:
            os.umask(old_umask)
        os.chmod(self.path, 0o600)
        sock.listen(8)
        self._sock = sock
        self._running = True
        self._thread = threading.Thread(target=self._accept_loop, name="easyhybrid-control", daemon=True)
        self._thread.start()
        events.set_active(self.dialogs_off)
        if self.exit_handlers:
            self._install_exit_handlers()
        try:
            with open(INFO_FILE, "w") as f:
                json.dump({"socket": self.path, "pid": os.getpid(),
                           "version": self.api._version(),
                           "started": time.strftime("%Y-%m-%d %H:%M:%S")}, f)
        except OSError:
            pass
        print("EasyHybrid control server listening on {}".format(self.path))
        self._sync_menu(True)
        return self.path

    def _sync_menu(self, active):
        """ Edit > AI Control Server follows the real state. """
        try:
            sync = getattr(self.main, "sync_control_menu", None)
            if sync is not None:
                sync(active)
        except Exception:
            pass

    @staticmethod
    def _remove_stale_sockets(folder):
        """ control-<pid>.sock of processes that no longer exist (killed
            instances cannot remove their own). """
        import re
        for name in os.listdir(folder):
            m = re.match(r"control(?:-(\d+))?\.sock$", name)
            if not m:
                continue
            path = os.path.join(folder, name)
            if m.group(1) is None:
                # . the default control.sock: stale when nothing answers
                if not _alive(path):
                    try:
                        os.unlink(path)
                    except OSError:
                        pass
                continue
            try:
                os.kill(int(m.group(1)), 0)
            except ProcessLookupError:
                try:
                    os.unlink(path)
                except OSError:
                    pass
            except OSError:
                pass

    def _install_exit_handlers(self):
        """ Removes the socket on normal exit and on SIGTERM (kill), which
            otherwise ends the process without running any cleanup. """
        import atexit
        import signal
        atexit.register(self.stop)
        try:
            import gi
            gi.require_version("Gtk", "3.0")
            from gi.repository import GLib, Gtk

            def on_sigterm():
                self.stop()
                Gtk.main_quit()
                return False

            GLib.unix_signal_add(GLib.PRIORITY_HIGH, signal.SIGTERM, on_sigterm)
        except Exception:
            pass

    def stop(self):
        if self._running:
            self._sync_menu(False)
        self._running = False
        events.set_active(False)
        try:
            if self._sock is not None:
                self._sock.close()
        finally:
            self._sock = None
            if os.path.exists(self.path):
                try:
                    os.unlink(self.path)
                except OSError:
                    pass
            try:
                with open(INFO_FILE) as f:
                    if json.load(f).get("pid") == os.getpid():
                        os.unlink(INFO_FILE)
            except (OSError, ValueError):
                pass

    @property
    def running(self):
        return self._running

    def apply_settings(self, settings):
        """ New AI settings while running (Edit > AI Control Settings...). The
            dict is updated in place: the API shares it. """
        self.settings.clear()
        self.settings.update(settings)
        self.dialogs_off = not self.settings.get("dialogs")
        if self._running:
            events.set_active(self.dialogs_off)

    # --------------------------------------------------------------- socket
    def _accept_loop(self):
        while self._running:
            try:
                conn, _ = self._sock.accept()
            except OSError:
                break
            threading.Thread(target=self._client_loop, args=(conn,), daemon=True).start()

    def _client_loop(self, conn):
        reader = conn.makefile("r", encoding="utf-8")
        try:
            for line in reader:
                line = line.strip()
                if not line:
                    continue
                try:
                    request = json.loads(line)
                except ValueError as error:
                    reply = {"id": None, "ok": False, "error": "Invalid JSON: {}".format(error)}
                else:
                    reply = self.handle(request)
                data = json.dumps(reply, default=str) + "\n"
                conn.sendall(data.encode("utf-8"))
        except OSError:
            pass
        finally:
            try:
                reader.close()
                conn.close()
            except OSError:
                pass

    # -------------------------------------------------------------- dispatch
    def handle(self, request):
        """ One request -> one reply (runs in a client thread). """
        rid = request.get("id")
        method = request.get("method")
        params = request.get("params") or {}
        timeout = float(request.get("timeout") or 60)
        seq0 = events.last_seq()
        reply = {"id": rid}
        try:
            if method == "describe":
                reply.update(ok=True, result=command_schemas())
            elif method == "wait":
                reply.update(ok=True, result=self._wait(**params))
            elif method == "result":
                reply.update(self._result(**params))
            elif method == "settings":
                reply.update(ok=True, result=dict(self.settings))
            else:
                if not method or method.startswith("_") or not callable(getattr(self.api, method, None)):
                    raise APIError("Unknown command '{}'. Use describe.".format(method))
                reply.update(self._call(method, params, timeout))
                wait = params.get("wait")
                if wait is None:
                    wait = self.settings.get("auto_wait")
                # . a remote job may take hours: never waited for in the call
                if method in SIMULATIONS and wait and reply.get("ok") and not reply.get("pending") \
                        and not (reply.get("result") or {}).get("remote_job"):
                    reply["result"] = self._finish(method, reply["result"])
        except APIError as error:
            reply.update(ok=False, error=str(error))
        except TypeError as error:
            reply.update(ok=False, error="Bad arguments for {}: {}".format(method, error))
        except Exception as error:
            reply.update(ok=False, error="{}: {}".format(type(error).__name__, error), traceback=traceback.format_exc())
        reply["events"] = events.since(seq0)
        if self.settings.get("compact_json"):
            reply = _strip_empty(reply)
            reply["format"] = "compact"
        return reply

    def _finish(self, method, job):
        """ auto-wait: waits for the job, then imports and summarizes its
            result (pes for scan/neb, us_check + wham for umbrella). """
        sys_id, job_id = job.get("system"), job.get("job")
        info = self._wait(sys=sys_id, job=job_id, timeout=self.settings.get("wait_limit", 600))
        out = {"job": job_id, "system": sys_id, "type": job.get("type"), "status": info.get("status")}
        if info.get("timeout"):
            out["note"] = "still running after {} s: call wait(sys={}, job={})".format(
                self.settings.get("wait_limit", 600), sys_id, job_id)
            return out
        if info.get("status") != "Finished" and not str(info.get("status", "")).startswith("Finished"):
            out["error"] = info.get("error")
            return out
        if info.get("warning"):
            out["warning"] = info["warning"]
        steps = []
        if method == "umbrella":
            steps = [("check", "us_check", {"sys": sys_id, "job": job_id}),
                     ("pmf", "wham", {"sys": sys_id, "job": job_id})]
        else:
            steps = [("imported", "import_result", {"sys": sys_id, "job": job_id})]
            if method in ("scan", "neb"):
                steps.append(("energies", "pes", {"sys": sys_id, "job": job_id}))
        for key, command, params in steps:
            r = self._call(command, params, 600)
            out[key] = r.get("result") if r.get("ok") else {"error": r.get("error")}
        return out

    def _run_in_main(self, fn, thread=False):
        """ Schedules fn in the GTK main thread (or, thread=True, in a
            worker thread: commands marked off_main_thread); returns the job dict. """
        job = {"done": threading.Event()}

        def run():
            try:
                job["result"] = fn()
            except APIError as error:
                job["error"] = str(error)
            except TypeError as error:
                job["error"] = "Bad arguments: {}".format(error)
                job["traceback"] = traceback.format_exc()
            except Exception as error:
                job["error"] = "{}: {}".format(type(error).__name__, error)
                job["traceback"] = traceback.format_exc()
            finally:
                job["done"].set()
            return False

        if thread:
            threading.Thread(target=run, name="easyhybrid-command", daemon=True).start()
        else:
            from gi.repository import GLib
            GLib.idle_add(run)
        return job

    def _main_value(self, fn, timeout=5.0):
        job = self._run_in_main(fn)
        if not job["done"].wait(timeout):
            return None
        return job.get("result")

    def _job_reply(self, job):
        if "error" in job:
            out = {"ok": False, "error": job["error"]}
            if job.get("traceback"):
                out["traceback"] = job["traceback"]
            return out
        return {"ok": True, "result": _jsonable(job.get("result"))}

    def _call(self, method, params, timeout):
        fn = getattr(self.api, method)
        job = self._run_in_main(lambda: fn(**params), thread=getattr(fn, "off_main_thread", False))
        start = time.time()
        while not job["done"].wait(0.5):
            elapsed = time.time() - start
            if elapsed > 1.0:
                dialogs = [d for d in (self._main_value(self.api.dialogs, 2.0) or []) if d.get("modal")]
                if dialogs:
                    return self._park(job, dialogs=dialogs,
                                      note="Waiting for a dialog: answer it with dialog_respond, "
                                           "then get the result with method 'result'.")
            if elapsed > timeout:
                return self._park(job, note="Still running after {:.0f} s; get the result with method 'result'.".format(timeout))
        return self._job_reply(job)

    def _park(self, job, **extra):
        rid = "r{}".format(next(self._rid))
        self._pending[rid] = job
        out = {"ok": True, "pending": rid}
        out.update(extra)
        return out

    def _result(self, request, timeout=30):
        job = self._pending.get(request)
        if job is None:
            raise APIError("No pending request '{}'.".format(request))
        if not job["done"].wait(float(timeout)):
            dialogs = [d for d in (self._main_value(self.api.dialogs, 2.0) or []) if d.get("modal")]
            return {"ok": True, "pending": request, "dialogs": dialogs}
        del self._pending[request]
        return self._job_reply(job)

    def _wait(self, sys=None, job=None, timeout=600, poll=1.0):
        """ Waits until the job leaves Queued/Running (or timeout). """
        deadline = time.time() + float(timeout)
        info = None
        while True:
            holder = self._run_in_main(lambda: self.api.job_status(sys=sys, job=job))
            if not holder["done"].wait(30):
                raise APIError("The interface is not responding.")
            if "error" in holder:
                raise APIError(holder["error"])
            info = holder["result"]
            if info.get("status") not in (None, "Queued", "Running", "Running..."):
                return info
            if time.time() > deadline:
                info["timeout"] = True
                return info
            time.sleep(float(poll))
