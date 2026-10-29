#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Event log of the automation layer. No GTK import: the dialog helpers
#  (simple_dialog.py) and the process manager call record() and active().
#
#  In automation mode (control server running with dialogs off) the
#  information/error dialogs and the "Finished!" dialog are not shown: their
#  text is recorded here and returned with the reply of the request that
#  produced it. Question dialogs (Yes/No) are still shown; the client sees
#  them with "dialogs" and answers with "dialog_respond".

import itertools
import threading
import time
from collections import deque

_lock    = threading.Lock()
_events  = deque(maxlen=500)
_counter = itertools.count(1)
_state   = {"active": False}


def set_active(flag):
    """ Turns automation mode on/off (dialogs recorded instead of shown). """
    _state["active"] = bool(flag)


def active():
    """ True while dialogs are recorded instead of shown. """
    return _state["active"]


def record(kind, text="", title=None, **data):
    """ Adds an event; returns its sequence number.

        kind: "info", "error", "question", "finished", "job_error", "log", ...
    """
    with _lock:
        seq = next(_counter)
        event = {"seq": seq, "time": time.strftime("%H:%M:%S"), "kind": kind,
                 "title": title, "text": str(text or "").strip()}
        event.update(data)
        _events.append(event)
    return seq


def last_seq():
    with _lock:
        return _events[-1]["seq"] if _events else 0


def since(seq=0, kinds=None):
    """ Events with sequence number > seq (optionally only some kinds). """
    with _lock:
        out = [dict(e) for e in _events if e["seq"] > seq]
    if kinds:
        kinds = set(kinds)
        out = [e for e in out if e["kind"] in kinds]
    return out
