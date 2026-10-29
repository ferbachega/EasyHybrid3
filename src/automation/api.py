#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  AutomationAPI (2026-10-28): the commands of the control server.
#
#  Every public method (no leading underscore) is one command: keyword
#  arguments in, a JSON-serialisable value out. They run in the GTK main
#  thread (the server schedules them there), so they may touch the
#  interface. The signature and the docstring are the documentation: the
#  server's "describe" command and the MCP server build their tool lists
#  from them, so keep the first docstring line short and document every
#  argument in an "Args:" block.
#
#  Conventions
#    sys    pDynamo system id (the "id" of `state`); None = active system
#    obj    object id (the "id" of an object in `state`); None = last
#           object of the system
#    atoms  0-based atom indices of the pDynamo system (as in the
#           reaction coordinate boxes): a list [0, 4, 7] or a string
#           "0-10,15"
#    frame  frame of the object; None = the last frame
#
#  Only standard-library imports at module level: the MCP server imports
#  this module (without GTK) to read the signatures.

import inspect
import io
import os
import re
import time
import contextlib

SIMULATION_TYPES = (
    "Energy", "Geometry_Optimization", "Relaxed_Surface_Scan", "Advanced_Relaxed_Surface_Scan",
    "Molecular_Dynamics", "Nudged_Elastic_Band", "Umbrella_Sampling", "Normal_Modes",
    "Energy_Refinement",
)

_RC_TYPES = {
    "distance": "simple_distance", "simple_distance": "simple_distance",
    "multiple_distance": "multiple_distance",
    "multiple_distance_4": "multiple_distance*4atoms", "multiple_distance*4atoms": "multiple_distance*4atoms",
    "dihedral": "dihedral",
}

_MNDO_METHODS = ("am1", "am1dphot", "pm3", "pm6", "rm1", "mndo")


class APIError(Exception):
    """ A clear, expected error (bad argument, missing object...): the
        server returns its message without a traceback. """


def _jsonable(value):
    """ numpy scalars/arrays, tuples, sets and objects -> JSON types. """
    try:
        import numpy as np
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
    except ImportError:
        pass
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _angle(a, b, c):
    """ Angle a-b-c in degrees. """
    import numpy as np
    u, v = a - b, c - b
    cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def _dihedral(a, b, c, d):
    """ Dihedral a-b-c-d in degrees, (-180, 180], pDynamo's sign. """
    import numpy as np
    b1, b2, b3 = b - a, c - b, d - c
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return -float(np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2))))


def off_main_thread(method):
    """ Marks a command that does not touch GTK and may take long (reading
        files, WHAM): the server runs it in its own thread, so the window
        stays responsive. """
    method.off_main_thread = True
    return method


def _prominence(F, i, sign=1.0):
    """ Prominence of the local maximum i of sign*F (minimum for sign=-1):
        how far it stands above the higher of the lowest points that separate
        it from a higher peak (or the end) on each side. """
    import numpy as np
    G = sign * np.asarray(F, dtype=float)
    lows = []
    for step in (-1, 1):
        j, low = i, G[i]
        while 0 <= j + step < len(G) and np.isfinite(G[j + step]) and G[j + step] <= G[i]:
            j += step
            low = min(low, G[j])
        lows.append(low)
    return float(G[i] - max(lows))


def _profile_1d(x, F, min_prominence=1.0):
    """ Summary of a 1D energy curve: extremes and interior local extrema
        that stand out by at least min_prominence (kJ/mol; smaller ones are
        noise). """
    import numpy as np
    x, F = np.asarray(x, dtype=float), np.asarray(F, dtype=float)
    ok = np.isfinite(F)
    idx = np.nonzero(ok)[0]
    if not len(idx):
        return {}
    imax, imin = int(idx[np.argmax(F[ok])]), int(idx[np.argmin(F[ok])])
    first, last = int(idx[0]), int(idx[-1])
    point = lambda i: {"point": int(i), "rc": round(float(x[i]), 4), "energy": round(float(F[i]), 3)}
    inner = [i for i in idx if first < i < last and np.isfinite(F[i - 1]) and np.isfinite(F[i + 1])]
    maxima = [dict(point(i), prominence=round(_prominence(F, i), 3)) for i in inner
              if F[i] > F[i - 1] and F[i] >= F[i + 1] and _prominence(F, i) >= min_prominence]
    minima = [dict(point(i), depth=round(_prominence(F, i, -1.0), 3)) for i in inner
              if F[i] < F[i - 1] and F[i] <= F[i + 1] and _prominence(F, i, -1.0) >= min_prominence]
    out = {"highest": point(imax), "lowest": point(imin), "first": point(first), "last": point(last),
           "first_to_last": round(float(F[last] - F[first]), 3),
           # . interior extrema: the barrier of a curve that ends uphill is not its highest point
           "local_maxima": maxima, "local_minima": minima}
    if maxima:
        top = max(maxima, key=lambda m: m["energy"])
        before = [F[i] for i in idx if i <= top["point"]]
        out["barrier_from_first"] = round(float(top["energy"] - F[first]), 3)
        out["barrier_from_lowest_before_top"] = round(float(top["energy"] - min(before)), 3)
    return out


def _drop_arrays(out, keys, keep_short=True, limit=25):
    """ Summary replies: data arrays and path cells left out (keep_short keeps
        1D arrays of up to `limit` values). """
    dropped = []
    for key in keys:
        value = out.get(key)
        if value is None:
            continue
        if keep_short and isinstance(value, list) and len(value) <= limit and not (value and isinstance(value[0], list)):
            continue
        out.pop(key)
        dropped.append(key)
    path = out.get("path")
    if isinstance(path, dict) and "cells" in path:
        path["cells"] = len(path["cells"])
    if dropped:
        out["omitted"] = "{} (detail='full' returns them)".format(", ".join(dropped))
    return out


def _surface_2d(F, rc1_of, rc2_of, path_from=None, path_to=None):
    """ Summary of a 2D surface F[i, j]: lowest/highest cells and the
        lowest-barrier path between two cells (default: the two minima). """
    import numpy as np
    from pdynamo.us_data import default_endpoints, lowest_barrier_path
    F = np.asarray(F, dtype=float)
    cell = lambda c: {"cell": [int(c[0]), int(c[1])], "rc1": round(float(rc1_of(*c)), 4),
                      "rc2": round(float(rc2_of(*c)), 4), "energy": round(float(F[c]), 3)}
    G = np.where(np.isfinite(F), F, np.inf)
    lo = np.unravel_index(int(np.argmin(G)), F.shape)
    hi = np.unravel_index(int(np.nanargmax(np.where(np.isfinite(F), F, -np.inf))), F.shape)
    out = {"lowest": cell(lo), "highest": cell(hi)}
    try:
        a, b = default_endpoints(F)
        if path_from is not None:
            a = tuple(int(v) for v in path_from)
        if path_to is not None:
            b = tuple(int(v) for v in path_to)
        path = lowest_barrier_path(F, a, b)
        if path:
            top = max(path, key=lambda c: F[c])
            out["path"] = {"from": cell(a), "to": cell(b), "top": cell(top),
                           "barrier": round(float(F[top] - F[a]), 3),
                           "reaction_energy": round(float(F[b] - F[a]), 3),
                           "cells": [[int(i), int(j)] for i, j in path],
                           "note": "path_from/path_to choose other end cells [i, j]"}
    except Exception as error:
        out["path_error"] = str(error)
    return out


def parse_indices(spec):
    """ [1, 2, 5] / "1-3,5" / 7 -> sorted list of ints. """
    if spec is None:
        return None
    if isinstance(spec, int):
        return [spec]
    if isinstance(spec, str):
        out = []
        for part in spec.replace(" ", "").split(","):
            if not part:
                continue
            if "-" in part[1:]:
                lo, hi = part.split("-", 1) if not part.startswith("-") else (part, None)
                out.extend(range(int(lo), int(hi) + 1))
            else:
                out.append(int(part))
        return out
    return [int(i) for i in spec]


class AutomationAPI:
    """ Commands of the EasyHybrid control server. """

    def __init__(self, main, settings=None):
        self.main = main
        self._command = None
        self.settings = settings          # . automation/settings.py; the server shares its own

    def _settings(self):
        if self.settings is None:
            from automation import settings as ai_settings
            self.settings = ai_settings.load()
        return self.settings

    def _summary(self, detail):
        """ True when the reply should be a summary (detail None = the AI settings). """
        if detail is None:
            return bool(self._settings().get("summaries"))
        return str(detail).lower() != "full"

    # ------------------------------------------------------------------ access
    @property
    def p_session(self):
        return self.main.p_session

    @property
    def vm_session(self):
        return self.main.vm_session

    @property
    def _cmd(self):
        """ The terminal's Command (show/hide/color/... and the 3.2 system
            commands), shared with the in-app terminal when it exists. """
        if self._command is None:
            from gui.windows.setup.easyhybrid_terminal import Command
            self._command = Command(self, self.vm_session)
        return self._command

    def _version(self):
        return getattr(self.main, "EASYHYBRID_VERSION", None) or getattr(self.main, "version", None)

    def _system_ids(self):
        """ Ids of the loaded systems (psystem keeps None placeholders). """
        return {k for k, v in self.p_session.psystem.items() if v is not None}

    def _system(self, sys=None):
        psystem = self.p_session.psystem
        if not self._system_ids():
            raise APIError("No system loaded. Use load or new_system first.")
        e_id = self.p_session.active_id if sys is None else int(sys)
        system = psystem.get(e_id)
        if system is None:
            raise APIError("No system with id {} (see state).".format(e_id))
        return e_id, system

    def _objects_of(self, e_id):
        return [v for v in self.vm_session.vm_objects_dic.values() if getattr(v, "e_id", None) == e_id]

    def _vobject(self, obj=None, sys=None):
        if obj is not None:
            vobject = self.vm_session.vm_objects_dic.get(int(obj))
            if vobject is None:
                raise APIError("No object with id {} (see state).".format(obj))
            return vobject
        e_id, _ = self._system(sys)
        objects = self._objects_of(e_id)
        if not objects:
            raise APIError("System {} has no object.".format(e_id))
        return objects[-1]

    def _frame_index(self, vobject, frame=None):
        n = len(vobject.frames)
        if n == 0:
            raise APIError("Object {} has no coordinates.".format(vobject.index))
        if frame is None or int(frame) == -1:
            return n - 1
        frame = int(frame)
        if not 0 <= frame < n:
            raise APIError("Object {} has frames 0-{}.".format(vobject.index, n - 1))
        return frame

    def _coords(self, vobject, frame=None):
        import numpy as np
        f = self._frame_index(vobject, frame)
        return np.asarray(vobject.frames[f], dtype=float).reshape(-1, 3)

    def _activate(self, e_id):
        """ Makes e_id the active system (radio button of the Objects panel). """
        if self.p_session.active_id == e_id:
            return
        for row in self.main.main_treeview.treestore:
            row[4] = (row[0] == e_id)
        self.p_session.active_id = e_id
        try:
            self.main.refresh_main_statusbar()
        except Exception:
            pass

    def _use_coordinates(self, vobject, frame=None):
        """ The pDynamo system takes the coordinates of vobject/frame. """
        self._activate(vobject.e_id)
        f = self._frame_index(vobject, frame)
        self.p_session.set_psystem_coordinates_from_vobject(vobject=vobject, system_id=vobject.e_id, frame=f)
        return f

    def _folder(self, folder, system):
        folder = folder or getattr(system, "e_working_folder", None) or os.getcwd()
        folder = os.path.abspath(os.path.expanduser(folder))
        if not os.path.isdir(folder):
            raise APIError("Folder does not exist: {}".format(folder))
        return folder

    def _unique_name(self, folder, name, suffix=".ptGeo"):
        """ name, or name_2, name_3... when folder/name+suffix exists. """
        candidate, n = name, 1
        while os.path.exists(os.path.join(folder, candidate + suffix)):
            n += 1
            candidate = "{}_{}".format(name, n)
        return candidate

    def _select_atoms(self, vobject, atoms=None, chain=None, resi=None, resn=None, name=None, symbol=None,
                      within=None, of=None, frame=None, by_residue=False):
        """ 0-based indices of the atoms of vobject that match. """
        import numpy as np
        if atoms is not None and str(atoms).strip().lower() == "all":
            atoms = None
            chosen = [a.index - 1 for a in vobject.atoms.values()]
        elif atoms is not None:
            chosen = parse_indices(atoms)
        elif any(v is not None for v in (chain, resi, resn, name, symbol)):
            chosen = [a.index - 1 for a in self._cmd._filter_atoms(vobject, chain=chain, resi=resi, resn=resn,
                                                                  name=name, symbol=symbol)]
        elif within is not None:
            chosen = [a.index - 1 for a in vobject.atoms.values()]
        else:
            chosen = []
        if within is not None:
            if of is None:
                raise APIError("within= needs of= (the reference atoms).")
            xyz = self._coords(vobject, frame)
            ref = xyz[parse_indices(of)]
            d2 = ((xyz[:, None, :] - ref[None, :, :]) ** 2).sum(-1).min(1)
            close = set(np.nonzero(d2 <= float(within) ** 2)[0].tolist())
            chosen = [i for i in chosen if i in close]
        if by_residue and chosen:
            residues = set()
            by_index = {a.index - 1: a for a in vobject.atoms.values()}
            for i in chosen:
                residue = getattr(by_index.get(i), "residue", None)
                if residue is not None:
                    residues.add(id(residue))
            chosen = [a.index - 1 for a in vobject.atoms.values() if id(getattr(a, "residue", None)) in residues]
        return sorted(set(chosen))

    def _job(self, system, job=None):
        history = system.e_job_history
        if not history:
            raise APIError("System {} has no job.".format(system.e_id))
        key = max(history) if job is None else int(job)
        if key not in history:
            raise APIError("System {} has no job {}.".format(system.e_id, job))
        return key, history[key]

    def _job_info(self, system, key, job):
        objects = [v.index for v in self._objects_of(system.e_id) if getattr(v, "key6", None) == job.get("key6")]
        error = job.get("error")
        info = {
            "system": system.e_id, "job": key, "type": job.get("simulation_type"),
            "status": job.get("status"), "potential": job.get("potential"),
            "logfile": job.get("logfile"), "warning": job.get("warning"),
            "error": None if error is None else "{}: {}".format(type(error).__name__, error),
            "object": objects[-1] if objects else None,
        }
        if job.get("remote"):
            info.update({"remote": job.get("remote"), "remote_job": job.get("remote_job"),
                         "remote_state": job.get("remote_state")})
        # . the process is still alive while 'status' is not set yet
        entry = getattr(self.p_session, "process_pool", {}).get(system.e_id)
        if info["status"] in (None, "Queued") and entry and entry[1].is_alive():
            info["status"] = "Running"
        return info

    @staticmethod
    def _remote_args(parameters, scope):
        """ remote / remote_cpus arguments of a simulation command -> parameters. """
        if scope.get("remote"):
            parameters["remote"] = scope["remote"]
        if scope.get("remote_cpus"):
            parameters["remote_cpus"] = int(scope["remote_cpus"])

    def _start(self, parameters, vobject, system):
        """ run_simulation() with the active system = vobject's system;
            returns the job record. """
        e_id = system.e_id
        from pdynamo.pDynamo2EasyHybrid.remote_mixin import REMOTE_TYPES
        remote = parameters.get("remote")
        if remote is None and parameters["simulation_type"] in REMOTE_TYPES:
            remote = getattr(self.p_session, "remote_target", None)
        if remote:
            parameters["remote"] = remote
            parameters.setdefault("obj1_key6", vobject.key6)
            step = system.e_step_counter
            record = self.p_session.run_simulation(parameters=parameters)
            if not record:
                raise APIError("The remote job was not prepared (see events).")
            return {"system": e_id, "job": step, "type": parameters["simulation_type"], "status": "Queued",
                    "remote": remote, "remote_job": record["id"], "results_folder": record["local_folder"],
                    "hint": "it runs on {}: wait(sys={}, job={}) or remote_jobs; results come back "
                            "automatically".format(remote, e_id, step)}
        entry = getattr(self.p_session, "process_pool", {}).get(e_id)
        if entry and entry[1].is_alive():
            raise APIError("System {} already has a job running (one job per system). Wait for it.".format(e_id))
        parameters.setdefault("obj1_key6", vobject.key6)
        step = system.e_step_counter
        # . run_simulation() returns False in every case: a new process_pool entry means started
        self.p_session.run_simulation(parameters=parameters)
        if getattr(self.p_session, "process_pool", {}).get(e_id) is entry:
            raise APIError("The job was not started (see events).")
        return {"system": e_id, "job": step, "type": parameters["simulation_type"],
                "status": "Running", "logfile": system.e_job_history.get(step, {}).get("logfile"),
                "hint": "wait(sys={}, job={}) then import_result(sys={}, job={})".format(e_id, step, e_id, step)}

    # ===================================================================== info
    def ping(self):
        """ Checks the connection; returns the EasyHybrid version. """
        return {"ok": True, "version": self._version(), "pid": os.getpid()}

    def state(self):
        """ Systems, objects, running jobs and open dialogs of the session.

        Start here: the ids returned are the sys= and obj= of the other
        commands.
        """
        from pdynamo import mm_inspect
        systems = []
        for e_id, s in self.p_session.psystem.items():
            if s is None:
                continue
            qc = None
            if getattr(s, "qcModel", None) is not None:
                state = getattr(s, "electronicState", None)
                qc = {"model": type(s.qcModel).__name__.replace("QCModel", ""),
                      "hamiltonian": getattr(s.qcModel, "hamiltonian", None),
                      "atoms": len(getattr(s, "e_qc_table", None) or []) or len(s.atoms),
                      "charge": getattr(state, "charge", None), "multiplicity": getattr(state, "multiplicity", None)}
            try:
                ff = mm_inspect.force_field_name(s)
            except Exception:
                ff = None
            systems.append({
                "id": e_id, "name": s.label, "tag": getattr(s, "e_tag", None), "atoms": len(s.atoms),
                "active": e_id == self.p_session.active_id, "force_field": ff, "qc": qc,
                "fixed_atoms": len(getattr(s, "e_fixed_table", None) or []),
                "restraints": [k for k, v in (getattr(s, "e_restraints_dict", None) or {}).items() if v and v[0]],
                "working_folder": getattr(s, "e_working_folder", None),
                "jobs": len(getattr(s, "e_job_history", None) or {}),
            })
        objects = [{"id": v.index, "name": v.name, "system": getattr(v, "e_id", None),
                    "frames": len(v.frames), "atoms": len(v.atoms), "visible": bool(v.active)}
                   for v in self.vm_session.vm_objects_dic.values()]
        running = []
        for e_id, entry in getattr(self.p_session, "process_pool", {}).items():
            if entry and entry[1].is_alive() and e_id in self.p_session.psystem:
                s = self.p_session.psystem[e_id]
                running.append({"system": e_id, "job": s.e_step_counter,
                                "type": s.e_job_history.get(s.e_step_counter, {}).get("simulation_type")})
        return {"version": self._version(), "systems": systems, "objects": objects,
                "frame": int(self.vm_session.get_frame() or 0), "running_jobs": running,
                "session_file": getattr(self.main, "session_filename", None), "dialogs": self.dialogs()}

    def atoms(self, obj=None, sys=None, atoms=None, chain=None, resi=None, resn=None, name=None, symbol=None,
              within=None, of=None, by_residue=False, frame=None, limit=200):
        """ Lists atoms (index, name, residue, chain, element, xyz, charge, QC/fixed).

        Args:
          obj, sys: object / system (defaults: last object of the active system)
          atoms: indices "0-10,15" or "all"
          chain, resi, resn, name, symbol: filters, e.g. resn="HIS", resi="10-20", symbol="C,N"
          within, of: keep atoms within this distance (A) of the atoms `of`
          by_residue: extend the selection to whole residues
          frame: coordinates of this frame (default last)
          limit: maximum number of atoms listed (the count is always given)
        """
        vobject = self._vobject(obj, sys)
        if atoms is None and all(v is None for v in (chain, resi, resn, name, symbol, within)):
            atoms = "all"
        chosen = self._select_atoms(vobject, atoms, chain, resi, resn, name, symbol, within, of, frame, by_residue)
        xyz = self._coords(vobject, frame)
        system = self.p_session.psystem.get(vobject.e_id)
        qc = set(getattr(system, "e_qc_table", None) or [])
        fixed = set(getattr(system, "e_fixed_table", None) or [])
        charges = None
        try:
            charges = system.mmState.charges
        except Exception:
            pass
        by_index = {a.index - 1: a for a in vobject.atoms.values()}
        rows = []
        for i in chosen[:int(limit)]:
            a = by_index.get(i)
            if a is None:
                continue
            rows.append({
                "index": i, "name": a.name, "element": getattr(a, "symbol", None),
                "resn": getattr(getattr(a, "residue", None), "name", None),
                "resi": getattr(getattr(a, "residue", None), "index", None),
                "chain": getattr(getattr(a, "chain", None), "name", None),
                "xyz": [round(float(c), 4) for c in xyz[i]],
                "charge": round(float(charges[i]), 4) if charges is not None else None,
                "qc": i in qc, "fixed": i in fixed,
            })
        return {"object": vobject.index, "system": vobject.e_id, "count": len(chosen),
                "indices": chosen if len(chosen) <= 5000 else None, "atoms": rows}

    def measure(self, atoms, obj=None, sys=None, frame=None, frames=None):
        """ Distance (2 atoms), angle (3) or dihedral (4), in A / degrees.

        Args:
          atoms: the 2-4 atom indices, or a list of such lists
          frame: one frame (default last)
          frames: "all" or "0-20" for a value per frame (trajectories, scans)
        """
        import numpy as np
        vobject = self._vobject(obj, sys)
        groups = atoms if atoms and isinstance(atoms[0], (list, tuple)) else [atoms]
        groups = [parse_indices(g) for g in groups]
        if frames is None:
            frame_list = [self._frame_index(vobject, frame)]
        elif str(frames).lower() == "all":
            frame_list = list(range(len(vobject.frames)))
        else:
            frame_list = parse_indices(frames)

        def value(xyz, g):
            p = [xyz[i] for i in g]
            if len(g) == 2:
                return float(np.linalg.norm(p[0] - p[1]))
            if len(g) == 3:
                return _angle(*p)
            if len(g) == 4:
                return _dihedral(*p)
            raise APIError("measure takes 2, 3 or 4 atoms, got {}.".format(len(g)))

        out = []
        for g in groups:
            kind = {2: "distance", 3: "angle", 4: "dihedral"}.get(len(g), "?")
            values = [round(value(self._coords(vobject, f), g), 4) for f in frame_list]
            out.append({"atoms": g, "kind": kind, "unit": "A" if kind == "distance" else "deg",
                        "values": values if frames is not None else None,
                        "value": values[0] if frames is None else None})
        return {"object": vobject.index, "frames": frame_list if frames is not None else frame_list[0], "results": out}

    # ================================================================ systems
    def load(self, path):
        """ Opens a file as File > Open (.easy session, .pkl/.yaml system, PDB/XYZ/MOL2...).

        Args:
          path: file path
        """
        path = os.path.abspath(os.path.expanduser(str(path)))
        if not os.path.isfile(path):
            raise APIError("File not found: {}".format(path))
        before_s = self._system_ids()
        before_o = set(self.vm_session.vm_objects_dic)
        self.vm_session.load(path)
        if path.endswith(".easy"):
            self.main.session_filename = path
        return {"new_systems": sorted(self._system_ids() - before_s),
                "new_objects": sorted(set(self.vm_session.vm_objects_dic) - before_o)}

    def new_system(self, coordinates, force_field="other", name=None, tag="molsys", topology=None,
                   parameters=None, parameter_folder=None):
        """ Creates a system as File > New (force field and coordinates).

        Args:
          coordinates: coordinate file (pdb, xyz, mol2, crd...) or a .pkl/.yaml system
          force_field: "amber" (topology=.prmtop), "charmm" (topology=.psf, parameters=[.prm...]),
                       "opls", "dyff", "pkl" (pDynamo file) or "other" (no force field)
          name, tag: system name and tag
          parameter_folder: OPLS/DYFF parameter set (default: the pDynamo OPLS protein set / dyff-1.0)
        """
        types = {"amber": 0, "charmm": 1, "opls": 2, "pkl": 3, "other": 4, "dyff": 5}
        if force_field not in types:
            raise APIError("force_field must be one of {}.".format(", ".join(types)))
        system_type = types[force_field]
        coordinates = os.path.abspath(os.path.expanduser(coordinates))
        if not os.path.isfile(coordinates):
            raise APIError("File not found: {}".format(coordinates))
        files = {"coordinates": coordinates, "charges": None}
        if system_type == 0:
            files["amber_prmtop"] = os.path.abspath(os.path.expanduser(topology or ""))
        elif system_type == 1:
            files["charmm_psf"] = os.path.abspath(os.path.expanduser(topology or ""))
            files["charmm_par"] = [os.path.abspath(os.path.expanduser(p)) for p in (parameters or [])]
        elif system_type == 2:
            files["prm_folder"] = parameter_folder or os.path.join(os.environ.get("PDYNAMO3_PARAMETERS", ""),
                                                                   "forceFields/opls/protein")
        elif system_type == 5:
            from pdynamo.dyff_edits import dyff_parameter_dir
            files["prm_folder"] = parameter_folder or dyff_parameter_dir("dyff-1.0")
            if coordinates.lower().endswith(".mol2"):
                try:
                    from util.file_parser import read_MOL2
                    files["charges"] = read_MOL2(coordinates)["charges"]
                except Exception:
                    files["charges"] = None
        before = self._system_ids()
        self.p_session.load_a_new_pDynamo_system_from_dict(input_files=files, system_type=system_type,
                                                           name=name, tag=tag)
        new = sorted(self._system_ids() - before)
        if not new:
            raise APIError("The system was not created (see events).")
        e_id = new[-1]
        return {"system": e_id, "name": self.p_session.psystem[e_id].label,
                "atoms": len(self.p_session.psystem[e_id].atoms),
                "objects": [v.index for v in self._objects_of(e_id)]}

    def save_session(self, path=None):
        """ Saves the session (.easy). Without path, to the current session file.

        Args:
          path: .easy file
        """
        path = path or getattr(self.main, "session_filename", None)
        if not path:
            raise APIError("No session file yet: give path=.../name.easy")
        path = os.path.abspath(os.path.expanduser(path))
        if not path.endswith(".easy"):
            path += ".easy"
        self.p_session.save_easyhybrid_session(filename=path)
        self.main.session_filename = path
        return {"saved": path}

    def activate(self, sys):
        """ Makes a system the active one (target of menus and of commands without sys=).

        Args:
          sys: system id
        """
        e_id, system = self._system(sys)
        self._activate(e_id)
        return {"active": e_id, "name": system.label}

    def set_qc(self, atoms=None, method="am1", charge=0, multiplicity=1, engine="pdynamo", sys=None, obj=None,
               chain=None, resi=None, resn=None, name=None, within=None, of=None, by_residue=False,
               restricted=True, orca_keywords="HF 6-31G*", gfn=2, cpus=1, scratch=None):
        """ Defines the QC region and the QC method (System > Set QC Region).

        Args:
          atoms: QC atoms ("0-15,20", a list, or "all" for a pure QC system); or use the filters below
          chain, resi, resn, name, within, of, by_residue: atom filters as in `atoms`
          method: pdynamo engine: am1, am1dphot, pm3, pm6, rm1, mndo
          engine: "pdynamo" (semiempirical), "orca", "xtb"
          charge, multiplicity: of the QC region
          orca_keywords: ORCA method line, e.g. "B3LYP 6-31G* D3"
          gfn: xTB GFN version (0, 1, 2)
          cpus: xTB threads
        """
        engine = str(engine).lower()
        if engine not in ("pdynamo", "orca", "xtb"):
            raise APIError("engine must be pdynamo, orca or xtb.")
        if engine == "pdynamo" and str(method).lower() not in _MNDO_METHODS:
            raise APIError("method must be one of {} (or engine=orca/xtb).".format(", ".join(_MNDO_METHODS)))
        vobject = self._vobject(obj, sys)
        e_id = vobject.e_id
        _, system = self._system(e_id)
        self._activate(e_id)
        if atoms is not None and str(atoms).strip().lower() == "all":
            qc_list = []
        else:
            qc_list = self._select_atoms(vobject, atoms, chain, resi, resn, name, None, within, of, None, by_residue)
            if not qc_list:
                raise APIError("No atom selected for the QC region.")
        if len(qc_list) == len(system.atoms):
            qc_list = []
        residue_dict = {}
        by_index = {a.index - 1: a for a in vobject.atoms.values()}
        for i in qc_list:
            residue = getattr(by_index.get(i), "residue", None)
            if residue is not None:
                residue_dict.setdefault(residue.index, []).append(i)
        system.e_qc_table = qc_list
        system.e_qc_residue_table = residue_dict

        parameters = {"charge": int(charge), "multiplicity": int(multiplicity), "isSpinRestricted": bool(restricted),
                      "qcmodel": "mndo", "method": str(method).lower(),
                      "energyTolerance": 0.0001, "densityTolerance": 1e-08, "maximumIterations": 550}
        scratch = scratch or os.environ.get("PDYNAMO3_SCRATCH") or os.path.expanduser("~")
        if engine == "pdynamo":
            if parameters["method"] not in _MNDO_METHODS:
                raise APIError("method must be one of {} (or engine=orca/xtb).".format(", ".join(_MNDO_METHODS)))
            parameters["qcengine"] = "pDynamo"
        elif engine == "orca":
            parameters.update({"qcengine": "ORCA", "orca_options": orca_keywords, "orca_scratch": scratch,
                               "random_scratch": False})
        elif engine == "xtb":
            parameters.update({"qcengine": "xTB", "gfn": int(gfn), "parallel": int(cpus), "acc": 1.0,
                               "iterations": 300, "fermi_temp": 300.0, "keywords": "",
                               "scratch": os.path.join(scratch, str(getattr(system, "e_tag", "") or "xtb")),
                               "xtb_backup_files": ["log"]})
        else:
            raise APIError("engine must be pdynamo, orca or xtb.")
        ok = self.p_session.define_a_new_QCModel(system=None, parameters=parameters, vismol_object=vobject)
        if not ok:
            raise APIError("The QC model was not defined (see events).")
        return {"system": e_id, "qc_atoms": len(qc_list) or len(system.atoms), "engine": parameters["qcengine"],
                "method": parameters["method"] if engine == "pdynamo" else (orca_keywords if engine == "orca" else "GFN{}-xTB".format(gfn)),
                "charge": int(charge), "multiplicity": int(multiplicity),
                "link_atoms": len(getattr(getattr(system, "qcState", None), "boundaryAtoms", None) or [])}

    def fix_atoms(self, atoms=None, free_within=None, of=None, by_residue=True, clear=False, sys=None, obj=None,
                  chain=None, resi=None, resn=None, name=None):
        """ Sets the fixed atoms of a system (System > Set as Fixed).

        Args:
          atoms: atoms to fix ("0-100", a list, or filters chain/resi/resn/name)
          free_within, of: instead, fix everything EXCEPT the atoms within this distance (A) of `of`
          by_residue: with free_within, keep whole residues free (default true)
          clear: free every atom
        """
        vobject = self._vobject(obj, sys)
        e_id = vobject.e_id
        _, system = self._system(e_id)
        self._activate(e_id)
        if clear:
            fixed = []
        elif free_within is not None:
            free = set(self._select_atoms(vobject, None, None, None, None, None, None, free_within, of, None, by_residue))
            fixed = [i for i in range(len(system.atoms)) if i not in free]
        else:
            fixed = self._select_atoms(vobject, atoms, chain, resi, resn, name)
            if not fixed:
                raise APIError("No atom selected.")
        self.p_session.define_free_or_fixed_atoms_from_iterable(fixed)
        for v in self._objects_of(e_id):
            self.p_session._apply_fixed_representation_to_vobject(system_id=None, vismol_object=v)
            self.p_session._apply_QC_representation_to_vobject(system_id=None, vismol_object=v)
        self.vm_session.vm_glcore.queue_draw()
        return {"system": e_id, "fixed_atoms": len(fixed), "free_atoms": len(system.atoms) - len(fixed)}

    def energy(self, sys=None, obj=None, frame=None):
        """ Single-point energy (kJ/mol) with its terms, at the coordinates of obj/frame.

        Args:
          obj, frame: coordinates to use (default: last frame of the last object)
        """
        vobject = self._vobject(obj, sys)
        _, system = self._system(vobject.e_id)
        if getattr(system, "mmModel", None) is None and getattr(system, "qcModel", None) is None:
            raise APIError("System {} has no energy model (force field or QC model).".format(system.e_id))
        f = self._use_coordinates(vobject, frame)
        energy = system.Energy(log=None)
        terms = getattr(getattr(system, "scratch", None), "energyTerms", None)
        return {"system": system.e_id, "object": vobject.index, "frame": f, "energy": round(float(energy), 6),
                "unit": "kJ/mol", "terms": {k: round(float(v), 6) for k, v in dict(terms or {}).items()}}

    # ============================================================= simulations
    def optimize(self, obj=None, sys=None, frame=None, optimizer="ConjugatedGradient", max_iterations=600,
                 rms_gradient=0.1, name=None, folder=None, save_trajectory=False, save_frequency=10, wait=None,
                 remote=None, remote_cpus=None):
        """ Geometry optimization (Simulate > Geometry Optimization). Returns the job; see wait.

        Args:
          obj, frame: starting coordinates
          optimizer: ConjugatedGradient, SteepestDescent, LFBGS, QuasiNewton, FIRE
          max_iterations, rms_gradient: convergence (rms gradient in kJ/mol/A)
          name: name of the result object / trajectory
          folder: output folder (default: the system's working folder)
          save_trajectory: also write the optimization path (.ptGeo)
          remote: run on this remote machine (see remotes); default: Simulate > Run Simulations On
          remote_cpus: CPUs requested from the remote machine (default: its profile)
          wait: wait for the job and return its result in this call (default: the AI settings)
        """
        vobject = self._vobject(obj, sys)
        _, system = self._system(vobject.e_id)
        self._use_coordinates(vobject, frame)
        folder = self._folder(folder, system)
        name = name or "{}_opt".format(system.e_tag or "system")
        parameters = {"simulation_type": "Geometry_Optimization", "trajectory_name": None, "dialog": False,
                      "folder": folder, "optimizer": optimizer, "maximumIterations": int(max_iterations),
                      "logFrequency": 10, "save_frequency": int(save_frequency),
                      "rmsGradientTolerance": float(rms_gradient), "save_format": None,
                      "save_traj": False, "save_pdb": False, "vobject_name": name}
        if save_trajectory:
            name = self._unique_name(folder, name)
            parameters.update({"save_traj": True, "trajectory_name": name + ".ptGeo", "save_format": ".ptGeo",
                               "vobject_name": name})
        self._remote_args(parameters, locals())
        return self._start(parameters, vobject, system)

    def _rc(self, rc, vobject, frame, system):
        """ Reaction coordinate dict of the scan/US windows from a short spec. """
        from gui.widgets.custom_widgets.reaction_coordinate import compute_sigma_a1_a3
        if not isinstance(rc, dict):
            raise APIError("A reaction coordinate is a dict: {type, atoms, step, nsteps, ...}.")
        rc_type = _RC_TYPES.get(str(rc.get("type", "distance")).lower())
        if rc_type is None:
            raise APIError("rc type must be distance, multiple_distance, multiple_distance_4 or dihedral.")
        atoms = parse_indices(rc.get("atoms"))
        need = {"simple_distance": 2, "multiple_distance": 3, "multiple_distance*4atoms": 4, "dihedral": 4}[rc_type]
        if not atoms or len(atoms) != need:
            raise APIError("{} needs {} atoms.".format(rc_type, need))
        xyz = self._coords(vobject, frame)
        import numpy as np
        d = lambda i, j: float(np.linalg.norm(xyz[i] - xyz[j]))
        s13, s31 = (1.0, -1.0) if rc_type.startswith("multiple") else (None, None)
        mass_weighted = bool(rc.get("mass_weighted", False)) and rc_type == "multiple_distance"
        if mass_weighted:
            s13, s31 = compute_sigma_a1_a3(vobject, atoms[0], atoms[2])
        if rc_type == "simple_distance":
            start = d(atoms[0], atoms[1])
        elif rc_type == "multiple_distance":
            start = s13 * d(atoms[0], atoms[1]) + s31 * d(atoms[1], atoms[2])
        elif rc_type == "multiple_distance*4atoms":
            start = d(atoms[0], atoms[1]) - d(atoms[2], atoms[3])
        else:
            start = _dihedral(*[xyz[i] for i in atoms])
        if rc.get("start") is not None:
            start = float(rc["start"])
        by_index = {a.index - 1: a for a in vobject.atoms.values()}
        out = {"rc_type": rc_type, "ATOMS": atoms, "ATOM_NAMES": [getattr(by_index.get(i), "name", "X") for i in atoms],
               "dminimum": round(start, 5), "sigma_pk1pk3": s13, "sigma_pk3pk1": s31,
               "nsteps": int(rc.get("nsteps", 10)), "dincre": float(rc.get("step", 0.1)),
               "force_constant": float(rc.get("force_constant", 4000.0))}
        if mass_weighted:
            out["MC"] = True
        return out

    def scan(self, rc1, rc2=None, obj=None, sys=None, frame=None, name="scan", folder=None, cpus=1,
             optimizer="ConjugatedGradient", max_iterations=1000, rms_gradient=0.1, wait=None,
             remote=None, remote_cpus=None):
        """ Relaxed surface scan, 1D or 2D (Simulate > Reaction Coordinate Scans).

        A reaction coordinate is a dict:
          {"type": "distance" | "multiple_distance" | "multiple_distance_4" | "dihedral",
           "atoms": [i, j(, k(, l))], "step": 0.1, "nsteps": 20,
           "force_constant": 4000, "mass_weighted": false, "start": null}
        multiple_distance (3 atoms) is d(i-j) - d(j-k), the atom j moving from i to k;
        multiple_distance_4 is d(i-j) - d(k-l). start defaults to the current value.

        Args:
          rc1: first coordinate
          rc2: second coordinate (2D scan, run on `cpus` processes)
          obj, frame: starting coordinates
          name: trajectory name (folder <name>.ptGeo)
          cpus: processes for a 2D scan
          remote: run on this remote machine (see remotes); default: Simulate > Run Simulations On
          remote_cpus: CPUs requested from the remote machine (default: its profile)
          wait: wait for the job and return its result in this call (default: the AI settings)
        """
        vobject = self._vobject(obj, sys)
        _, system = self._system(vobject.e_id)
        f = self._use_coordinates(vobject, frame)
        folder = self._folder(folder, system)
        name = self._unique_name(folder, name)
        parameters = {"simulation_type": "Relaxed_Surface_Scan", "logFrequency": 50, "optimizer": optimizer,
                      "folder": folder, "maximumIterations": float(max_iterations),
                      "rmsGradientTolerance": float(rms_gradient), "traj_folder_name": name, "vobject_name": name,
                      "initial_coordinates": vobject.name, "_is_ts_centered": False,
                      "RC1": self._rc(rc1, vobject, f, system)}
        if rc2:
            parameters.update({"RC2": self._rc(rc2, vobject, f, system), "NmaxThreads": int(cpus),
                               "traj_type": "pklfolder2D"})
        else:
            parameters.update({"RC2": None, "NmaxThreads": 1, "traj_type": "pklfolder"})
        self._remote_args(parameters, locals())
        job = self._start(parameters, vobject, system)
        job["rc1_start"] = parameters["RC1"]["dminimum"]
        if rc2:
            job["rc2_start"] = parameters["RC2"]["dminimum"]
        return job

    def md(self, obj=None, sys=None, frame=None, steps=1000, time_step=0.001, integrator="Langevin",
           temperature=300.0, collision_frequency=25.0, name=None, folder=None, save_frequency=50,
           log_frequency=100, seed=None, temperature_end=None, scale_option="constant", scale_frequency=10,
           temperature_coupling=0.1, pressure=None, pressure_coupling=2000.0, wait=None,
           remote=None, remote_cpus=None):
        """ Molecular dynamics (Simulate > Molecular Dynamics). Returns the job; see wait.

        Args:
          steps, time_step: number of steps and step in ps (0.001 = 1 fs)
          integrator: Langevin, Verlet or LeapFrog
          temperature: K; with temperature_end and scale_option "linear"/"exponential" a heating ramp
          collision_frequency: Langevin, ps^-1
          name: trajectory name (None = no trajectory saved)
          save_frequency: steps between saved frames
          pressure: atm, LeapFrog only (NPT)
          remote: run on this remote machine (see remotes); default: Simulate > Run Simulations On
          remote_cpus: CPUs requested from the remote machine (default: its profile)
          wait: wait for the job and return its result in this call (default: the AI settings)
        """
        import random
        vobject = self._vobject(obj, sys)
        _, system = self._system(vobject.e_id)
        self._use_coordinates(vobject, frame)
        folder = self._folder(folder, system)
        integrator = {"langevin": "Langevin", "verlet": "Verlet", "leapfrog": "LeapFrog"}.get(str(integrator).lower())
        if integrator is None:
            raise APIError("integrator must be Langevin, Verlet or LeapFrog.")
        leapfrog = integrator == "LeapFrog"
        parameters = {
            "simulation_type": "Molecular_Dynamics", "integrator": integrator, "logFrequency": int(log_frequency),
            "seed": int(seed) if seed is not None else random.randint(1, 10 ** 6), "steps": int(steps),
            "timeStep": float(time_step), "trajectories": None, "folder": folder,
            "temperatureScaleFrequency": int(scale_frequency), "temperatureScaleOption": scale_option,
            "temperatureStart": float(temperature),
            "temperatureStop": None if scale_option == "constant" else float(temperature_end or temperature),
            "pressure": float(pressure) if (leapfrog and pressure is not None) else None,
            "pressureControl": bool(leapfrog and pressure is not None),
            "pressureCoupling": float(pressure_coupling) if leapfrog else None,
            "temperature": float(temperature), "temperatureControl": True,
            "temperatureCoupling": float(temperature_coupling) if leapfrog else None,
            "collisionFrequency": float(collision_frequency) if integrator == "Langevin" else None,
            "restraints": None, "trajectory_name": None,
        }
        if name:
            name = self._unique_name(folder, name)
            parameters.update({"trajectory_name": name, "trajectory_format": ".ptGeo",
                               "trajectory_frequency": int(save_frequency)})
        self._remote_args(parameters, locals())
        return self._start(parameters, vobject, system)

    def neb(self, reactant, product, reactant_frame=None, product_frame=None, images=11, max_iterations=600,
            rms_gradient=0.1, spring=500.0, spline=True, name="neb", folder=None, cpus=1, wait=None,
            remote=None, remote_cpus=None):
        """ Nudged elastic band between two objects of the same system (Simulate > NEB).

        Args:
          reactant, product: object ids (the last frame of each unless *_frame is given)
          images: number of structures, ends included
          spring: spring force constant
          spline: spline redistribution of the images
          cpus: processes for the image energies
          remote: run on this remote machine (see remotes); default: Simulate > Run Simulations On
          remote_cpus: CPUs requested from the remote machine (default: its profile)
          wait: wait for the job and return its result in this call (default: the AI settings)
        """
        import copy
        v1, v2 = self._vobject(reactant), self._vobject(product)
        if v1.e_id != v2.e_id:
            raise APIError("Reactant and product must belong to the same system.")
        _, system = self._system(v1.e_id)
        folder = self._folder(folder, system)
        name = self._unique_name(folder, name)
        self._use_coordinates(v2, product_frame)
        prod = copy.deepcopy(system.coordinates3)
        self._use_coordinates(v1, reactant_frame)
        reac = copy.deepcopy(system.coordinates3)
        parameters = {
            "simulation_type": "Nudged_Elastic_Band", "number_of_structures": int(images), "logFrequency": 1,
            "maximumIterations": int(max_iterations), "rmsGradientTolerance": float(rms_gradient),
            "trajectory_name": name, "folder": folder, "vobject_name": name,
            "fixedTerminalImages": True, "forceOneSingleImageOptimization": False,
            "forceSingleImageOptimizations": False, "forceSplineRedistributionCheckPerIteration": False,
            "freezeRMSGradientTolerance": 0.0, "optimizer": None, "poolFactory": int(cpus),
            "rmsGradientToleranceScale": 0.25, "splineRedistributionTolerance": 1.5,
            "springForceConstant": float(spring), "useSplineRedistribution": bool(spline),
            "reac_coordinates": reac, "prod_coordinates": prod, "obj1_key6": v1.key6, "obj2_key6": v2.key6,
            "external_coords": False,
        }
        self._remote_args(parameters, locals())
        return self._start(parameters, v1, system)

    def simulate(self, parameters, obj=None, sys=None, frame=None):
        """ Any simulation from its full parameter dictionary (advanced).

        The dictionary is the one the Simulate windows build (see the
        "Export Script" output or pdynamo/p_methods). "simulation_type" is
        required; "folder" defaults to the working folder.

        Args:
          parameters: the dictionary
          obj, frame: starting coordinates
        """
        if not isinstance(parameters, dict) or "simulation_type" not in parameters:
            raise APIError("parameters must be a dict with simulation_type.")
        vobject = self._vobject(obj, sys)
        _, system = self._system(vobject.e_id)
        self._use_coordinates(vobject, frame)
        parameters = dict(parameters)
        parameters["folder"] = self._folder(parameters.get("folder"), system)
        self._remote_args(parameters, locals())
        return self._start(parameters, vobject, system)

    # ==================================================================== jobs
    def jobs(self, sys=None):
        """ Jobs of a system (or of all systems with sys="all"), newest last.

        Args:
          sys: system id, or "all"
        """
        if str(sys).lower() == "all":
            systems = [s for s in self.p_session.psystem.values() if s is not None]
        else:
            systems = [self._system(sys)[1]]
        out = []
        for s in systems:
            for key in sorted(s.e_job_history):
                out.append(self._job_info(s, key, s.e_job_history[key]))
        return {"jobs": out}

    def job_status(self, sys=None, job=None):
        """ Status of one job (default: the last job of the system).

        Args:
          sys, job: system id and job number
        """
        _, system = self._system(sys)
        key, job_record = self._job(system, job)
        return self._job_info(system, key, job_record)

    def import_result(self, sys=None, job=None, name=None):
        """ Imports the trajectory/frames of a finished job as a new object (the "Import..." of Finished!).

        Single-structure results (optimization without trajectory) already
        are an object: their id is returned.

        Args:
          sys, job: system id and job number (default: last job)
          name: name of the new object
        """
        _, system = self._system(sys)
        key, job_record = self._job(system, job)
        info = self._job_info(system, key, job_record)
        if info["status"] not in ("Finished", "Finished (not converged)"):
            raise APIError("Job {} is {}.".format(key, info["status"]))
        if job_record.get("simulation_type") == "Umbrella_Sampling":
            raise APIError("An umbrella sampling is analyzed with us_check and wham (sys={}, job={}).".format(
                system.e_id, key))
        logfile = job_record.get("logfile")
        folder = os.path.dirname(logfile) if logfile else None
        if not folder or not folder.endswith(".ptGeo") or not os.path.isdir(folder):
            if info["object"] is not None:
                return {"object": info["object"], "frames": len(self.vm_session.vm_objects_dic[info["object"]].frames),
                        "note": "single structure result"}
            raise APIError("Job {} wrote no trajectory to import.".format(key))
        data_type = "pklfolder"
        if any(re.match(r"frame\d+_\d+\.pkl$", f) for f in os.listdir(folder)):
            data_type = "pklfolder2D"
        try:
            from pdynamo.LogFileWriter import LogFileReader
            if LogFileReader(logfile).type is None:
                logfile = None
        except Exception:
            logfile = None
        name = name or "{}_{}_{}".format(system.label.strip(), job_record.get("simulation_type"), key)
        parameters = {"system_id": system.e_id, "data_path": folder, "data_type": data_type,
                      "new_vobj_name": name, "vobject_id": None, "vobject": None, "logfile": logfile,
                      "isAppend": False, "first": None, "last": None, "stride": None}
        self.p_session.import_data(parameters)
        vobject = parameters.get("vobject")
        if vobject is None:
            raise APIError("Import failed (see events).")
        return {"object": vobject.index, "name": vobject.name, "frames": len(vobject.frames),
                "data_type": data_type, "energies": bool(logfile)}

    def pes(self, sys=None, job=None, logfile=None, relative=True, path_from=None, path_to=None, detail=None):
        """ Energy profile/surface of a scan, NEB or refinement log (numbers, not a plot).

        1D: RC values, energies, extremes, interior local maxima (barriers)
        and minima. 2D: matrices with cell [i, j] = RC1 step i, RC2 step j,
        lowest/highest cells and the lowest-barrier path between the two
        minima (check it: isolated artefact cells can attract it; then give
        path_from/path_to).

        Args:
          sys, job: job whose log to read (default: last job of the active system)
          logfile: or the path of an output.log
          relative: energies relative to the lowest point (default true)
          path_from, path_to: 2D only, end cells [i, j] of the path (default: the two minima)
          detail: "summary" (no data arrays) or "full" (default: the AI settings)
        """
        from pdynamo.LogFileWriter import LogFileReader
        if logfile is None:
            _, system = self._system(sys)
            _, job_record = self._job(system, job)
            logfile = job_record.get("logfile")
        if not logfile or not os.path.isfile(logfile):
            raise APIError("Log file not found: {}".format(logfile))
        data = LogFileReader(logfile).get_data()
        if not data:
            raise APIError("No energy table in {}".format(logfile))
        import numpy as np
        Z = np.asarray(data["Z"], dtype=float)
        if relative:
            Z = Z - np.nanmin(Z)
        out = {"logfile": logfile, "type": data.get("type"), "unit": "kJ/mol",
               "rc1_atoms": data.get("RC1_indexes"), "rc2_atoms": data.get("RC2_indexes")}
        if data.get("type") == "plot2D":
            # . LogFileReader stores Z[j][i] (rows = RC2): transposed so that cell [i, j] is
            #   (RC1 step i, RC2 step j), the frame{i}_{j}.pkl of the scan
            Z = Z.T
            rc1, rc2 = np.asarray(data["RC1"], dtype=float).T, np.asarray(data["RC2"], dtype=float).T
            out.update({"shape": list(Z.shape), "RC1": np.round(rc1, 4).tolist(), "RC2": np.round(rc2, 4).tolist(),
                        "Z": np.round(Z, 3).tolist(),
                        "note": "cell [i, j] = RC1 step i, RC2 step j (structure frame{i}_{j}.pkl of the scan)"})
            out.update(_surface_2d(Z, lambda i, j: rc1[i, j], lambda i, j: rc2[i, j], path_from, path_to))
        else:
            rc1 = np.asarray(data["RC1"], dtype=float)
            out.update({"points": len(Z), "RC1": np.round(rc1, 4).tolist(), "Z": np.round(Z, 3).tolist()})
            out.update(_profile_1d(rc1, Z))
        if self._summary(detail):
            _drop_arrays(out, ("RC1", "RC2", "Z"), keep_short=out.get("type") != "plot2D")
        return out

    @off_main_thread
    def job_summary(self, sys=None, job=None, logfile=None, detail=None):
        """ What the log of a job says, as numbers: the "Plots" tab of View Log.

        Optimization: status, iterations, final energy and RMS gradient. MD:
        means of T and energies, drift. Scans: barriers, minima, 2D path.
        NEB: highest image, barrier. Umbrella sampling: overlap/deviation
        verdict. Normal modes: imaginary frequencies. Plus the ORCA/xTB files
        next to the log (energies of each calculation, last SCF, gap).

        Args:
          sys, job: the job (default: the last job of the active system)
          logfile: or the path of a log
          detail: "summary" (no plot data) or "full" (default: the AI settings)
        """
        from util.log_analysis import analyze
        sim_type = None
        if logfile is None:
            _, system = self._system(sys)
            _, job_record = self._job(system, job)
            logfile = job_record.get("logfile")
            sim_type = job_record.get("simulation_type")
        if not logfile or not os.path.isfile(logfile):
            raise APIError("Log file not found: {}".format(logfile))
        result = analyze(logfile, sim_type)
        out = {"logfile": logfile, "kind": result["kind"], "summary": dict(result["summary"])}
        if result.get("qc"):
            out["qc"] = {"program": result["qc"]["program"], "summary": dict(result["qc"]["summary"])}
        if self._summary(detail):
            out["plots"] = [p["title"] for p in result["panels"]]
        else:
            out["plots"] = result["panels"]
            if result.get("qc"):
                out["qc"]["plots"] = result["qc"]["panels"]
        return out

    # ========================================================== umbrella sampling
    def _us_folder(self, folder=None, sys=None, job=None):
        """ Folder of an umbrella sampling run: given, or that of a job. """
        if folder is None:
            _, system = self._system(sys)
            _, job_record = self._job(system, job)
            if job_record.get("simulation_type") != "Umbrella_Sampling":
                raise APIError("Job {} is a {}, not an umbrella sampling: give folder=.".format(
                    job_record.get("step_counter", job), job_record.get("simulation_type")))
            folder = os.path.dirname(job_record.get("logfile") or "")
        folder = os.path.abspath(os.path.expanduser(str(folder)))
        if not os.path.isdir(os.path.join(folder, "data_collection")):
            raise APIError("{} is not an umbrella sampling folder (no data_collection/).".format(folder))
        return folder

    def umbrella(self, source=None, sys=None, job=None, obj=None, frame=None, rc1=None, rc2=None,
                 force_constant=None, cpus=1, steps_eq=20000, steps_dc=50000, time_step=0.001,
                 integrator="Langevin", temperature=300.0, collision_frequency=25.0, temperature_coupling=0.1,
                 optimize=False, opt_max_iterations=600, opt_rms_gradient=0.1, save_frequency=100,
                 log_frequency=10, name="umbrella", folder=None, seed=None, wait=None,
                 remote=None, remote_cpus=None):
        """ Umbrella sampling (Simulate > Umbrella Sampling). Returns the job; see wait, us_check, wham.

        Two ways to define the windows:
          * from a scan trajectory (recommended; windows run in parallel):
            source=<folder .ptGeo> or sys/job of a finished scan. Each frame
            is a window, its target the RC value of the frame; the reaction
            coordinate is read from the scan's output.log (rc1/rc2 override it).
          * from one structure: obj (+frame) and rc1 with step and nsteps,
            the windows made one after the other (serial).
        Spacing and force constant go together: windows d apart overlap when
        d is about 2-3 sigma, sigma = sqrt(kT/2k) (k=600: 0.046 A at 300 K).

        Args:
          source: scan trajectory folder (.ptGeo); or sys and job of the scan
          obj, frame: starting structure (mode "from one structure")
          rc1, rc2: reaction coordinates as in scan; optional from a trajectory
          force_constant: restraint force constant, kJ/mol/A^2 (or per rc: "force_constant")
          cpus: windows run at the same time (from a trajectory)
          steps_eq, steps_dc: equilibration and data-collection steps per window
          time_step: ps; integrator: Langevin, Verlet, LeapFrog; temperature: K
          optimize: optimize each window with its restraint first
          save_frequency: steps between saved frames
          name: folder of the run (inside folder)
          remote: run on this remote machine (see remotes); default: Simulate > Run Simulations On
          remote_cpus: CPUs requested from the remote machine (default: its profile)
          wait: wait for the job and return its result in this call (default: the AI settings)
        """
        import random
        md_integrator = {"langevin": "Langevin", "verlet": "Verlet", "leapfrog": "LeapFrog"}.get(str(integrator).lower())
        if md_integrator is None:
            raise APIError("integrator must be Langevin, Verlet or LeapFrog.")
        from_trajectory = source is not None or job is not None
        if from_trajectory:
            if source is None:
                _, system = self._system(sys)
                _, job_record = self._job(system, job)
                source = os.path.dirname(job_record.get("logfile") or "")
            source = os.path.abspath(os.path.expanduser(str(source)))
            if not os.path.isdir(source):
                raise APIError("Trajectory folder not found: {}".format(source))
            frames = sorted(f for f in os.listdir(source) if re.match(r"frame\d+(_\d+)?\.pkl$", f))
            if not frames:
                raise APIError("No frame*.pkl in {} (a scan .ptGeo folder is expected).".format(source))
            vobject = self._vobject(obj, sys)
        else:
            if rc1 is None:
                raise APIError("Give source= (a scan trajectory) or obj= and rc1= with step and nsteps.")
            vobject = self._vobject(obj, sys)
            frames = []
        _, system = self._system(vobject.e_id)
        f = self._use_coordinates(vobject, frame)
        folder = self._folder(folder, system)
        name = self._unique_name(folder, name, suffix="")

        def with_k(rc_dict, spec):
            k = (spec or {}).get("force_constant") if isinstance(spec, dict) else None
            k = k if k is not None else force_constant
            if k is None:
                raise APIError("Give force_constant= (kJ/mol/A^2 or kJ/mol/rad^2 for a dihedral).")
            rc_dict["force_constant"] = float(k)
            return rc_dict

        detected = {}
        if from_trajectory and (rc1 is None or (rc2 is None and any("_" in fr[5:] for fr in frames))):
            from pdynamo.LogFileWriter import detect_reaction_coordinates_from_log
            detected = detect_reaction_coordinates_from_log(os.path.join(source, "output.log")) or {}

        def coordinate(spec, key):
            if spec is not None:
                return with_k(self._rc(spec, vobject, f, system), spec)
            data = detected.get(key)
            if not data:
                return None
            if data.get("rc_type") == "advanced":
                return with_k({"rc_type": "advanced", "RC": [list(r) for r in data["RC"]],
                               "nsteps": 0, "dincre": 0.0, "dminimum": 0.0}, None)
            out = {"rc_type": data["rc_type"], "ATOMS": list(data["ATOMS"]),
                   "ATOM_NAMES": list(data.get("ATOM_NAMES") or []),
                   "sigma_pk1pk3": data.get("sigma_pk1pk3", 1.0 if data["rc_type"].startswith("multiple") else None),
                   "sigma_pk3pk1": data.get("sigma_pk3pk1", -1.0 if data["rc_type"].startswith("multiple") else None),
                   "dminimum": data.get("dminimum") or 0.0, "nsteps": 0, "dincre": 0.0}
            return with_k(out, None)

        RC1 = coordinate(rc1, "RC1")
        if RC1 is None:
            raise APIError("No reaction coordinate found in {}/output.log: give rc1=.".format(source))
        RC2 = coordinate(rc2, "RC2") if (rc2 is not None or detected.get("RC2")) else None
        two_d = RC2 is not None
        if from_trajectory and two_d != any(re.match(r"frame\d+_\d+\.pkl$", fr) for fr in frames):
            raise APIError("The trajectory is {} but {} reaction coordinate(s) were given.".format(
                "2D" if not two_d else "1D", 2 if two_d else 1))
        md = {"simulation_type": "Umbrella_Sampling", "integrator": md_integrator, "logFrequency": int(log_frequency),
              "seed": int(seed) if seed is not None else random.randint(1, 10 ** 6), "normal_deviate_generator": None,
              "steps_eq": int(steps_eq), "steps_dc": int(steps_dc), "timeStep": float(time_step),
              "trajectory_frequency": int(save_frequency), "trajectory_frequency_dc_ptRes": 1,
              "trajectory_frequency_dc_ptGeo": int(save_frequency), "temperatureScaleFrequency": 10,
              "temperatureStart": float(temperature), "pressure": None, "pressureControl": False,
              "pressureCoupling": None, "temperature": float(temperature), "temperatureControl": True,
              "temperatureCoupling": float(temperature_coupling) if md_integrator != "Verlet" else None,
              "collisionFrequency": float(collision_frequency) if md_integrator == "Langevin" else None}
        parameters = {
            "simulation_type": "Umbrella_Sampling", "input_type": 1 if from_trajectory else 0,
            "vobj": None if from_trajectory else vobject.name, "source_folder": source if from_trajectory else None,
            "folder": folder, "RC1": RC1, "RC2": RC2,
            "NmaxThreads": int(cpus) if (from_trajectory or two_d) else 1,
            "OPT_parm": {"optimizer": "ConjugatedGradient", "maximumIterations": int(opt_max_iterations),
                         "rmsGradientTolerance": float(opt_rms_gradient)} if optimize else None,
            "MD_parm": md, "traj_folder_name": name,
        }
        self._remote_args(parameters, locals())
        job_info = self._start(parameters, vobject, system)
        job_info["windows"] = len(frames) if from_trajectory else None
        job_info["us_folder"] = os.path.join(folder, name)
        job_info["hint"] = "wait(sys={0}, job={1}) then us_check(sys={0}, job={1}) and wham(sys={0}, job={1})".format(
            system.e_id, job_info["job"])
        return job_info

    @off_main_thread
    def us_check(self, folder=None, sys=None, job=None, halves=True, detail=None):
        """ Diagnostics of an umbrella sampling run, as the Umbrella Sampling Browser computes them.

        Per window: target, mean, std, deviation (|mean - target| in units of
        the restraint width sqrt(kT/2k); above 3 the window slid away), drift
        between halves, samples. Neighbour overlaps (below 0.03 WHAM cannot
        connect them). With halves, the PMF of each half of the data
        (own WHAM) and their largest difference (convergence).

        Args:
          folder: the run folder (or sys and job of the umbrella sampling)
          halves: compute the convergence test (a few seconds 1D, about a minute 2D)
          detail: "summary" (only the windows and pairs with problems) or "full" (default: the AI settings)
        """
        import numpy as np
        from pdynamo.us_data import USRun
        folder = self._us_folder(folder, sys, job)
        run = USRun(folder)
        stats = run.stats()
        windows, deviating = [], []
        for key in run.keys:
            st = stats.get(key) or {}
            row = {"window": list(key), "samples": st.get("n", 0)}
            if st.get("n"):
                row.update({"target": None if st.get("target") is None else np.round(st["target"], 4).tolist(),
                            "mean": np.round(st["mean"], 4).tolist(), "std": np.round(st["std"], 4).tolist(),
                            "drift": np.round(st["drift"], 4).tolist(),
                            "deviation": round(float(st["deviation"]), 2)})
                if st["deviation"] > 3.0:
                    deviating.append(list(key))
            windows.append(row)
        pairs = run.overlap_1d()[0] if run.dims == 1 else run.overlap_2d()
        overlaps = [{"pair": [list(a), list(b)], "overlap": round(o, 3)} for a, b, o in pairs]
        low = [p for p in overlaps if p["overlap"] < 0.03]
        out = {"folder": folder, "dims": run.dims, "shape": list(run.shape), "windows_found": len(run.keys),
               "temperature": run.temperature, "windows": windows,
               "overlaps": overlaps, "min_overlap": min((p["overlap"] for p in overlaps), default=None),
               "low_overlap_pairs": low, "deviating_windows": deviating}
        issues, warnings = [], []
        if deviating:
            issues.append("{} window(s) more than 3 sigma from their target: the restraint did not hold the "
                          "coordinate (raise the force constant, or add windows on the steep part)".format(len(deviating)))
        if low:
            issues.append("{} neighbour pair(s) with overlap < 0.03: gaps WHAM cannot bridge (add windows or "
                          "lower the force constant)".format(len(low)))
        if halves:
            centers, F1 = run.wham(part=(0.0, 0.5))
            _, F2 = run.wham(part=(0.5, 1.0))
            D = np.abs(np.asarray(F2) - np.asarray(F1))
            whole = float(np.nanmax(D)) if np.isfinite(D).any() else None
            # . the histogram ends are thinly sampled: judge between the first and last targets
            largest = whole
            if run.dims == 1:
                targets = sorted(t[0] for t in (run.target(k) for k in run.keys) if t)
                inside = (np.asarray(centers) >= targets[0]) & (np.asarray(centers) <= targets[-1]) & np.isfinite(D)
                if inside.any():
                    largest = float(D[inside].max())
            out["halves_largest_difference"] = None if whole is None else round(whole, 2)
            out["halves_difference_within_targets"] = None if largest is None else round(largest, 2)
            if largest is not None and largest > 4.0:
                warnings.append("the PMFs of the two halves differ by up to {:.1f} kJ/mol: compare it with the "
                                "barrier; a longer data collection reduces it".format(largest))
        out["issues"] = issues
        out["warnings"] = warnings
        out["verdict"] = "fail" if issues else ("warning" if warnings else "OK")
        if self._summary(detail):
            bad = {tuple(k) for k in deviating}
            out["windows"] = [w for w in windows if tuple(w["window"]) in bad]
            out["overlaps"] = low
            out["note"] = "summary: only windows and pairs with problems; detail='full' lists all"
        return out

    @off_main_thread
    def wham(self, folder=None, sys=None, job=None, bins=100, bins_rc2=None, temperature=None,
             max_iterations=1000, rms_gradient=0.001, resamples=0, name="wham", path_from=None, path_to=None,
             detail=None):
        """ PMF of an umbrella sampling run with pDynamo's WHAM (Analysis > WHAM), as numbers.

        Writes <run>/wham/<name>_pmf.log. 1D: RC, PMF, extremes, local
        maxima (barriers) and minima. 2D: matrix (rows = RC1, columns = RC2)
        and the lowest-barrier path between the two minima.

        Args:
          folder: the run folder (or sys and job of the umbrella sampling)
          bins: histogram bins (RC1); bins_rc2: RC2 bins of a 2D run (default = bins)
          temperature: K (default: that of the run)
          resamples: bootstrap resamples for error bars (0 = none)
          path_from, path_to: 2D only, end cells [i, j] of the path
          detail: "summary" (no PMF arrays; they are in pmf_file) or "full" (default: the AI settings)
        """
        import glob
        import numpy as np
        from pdynamo.p_methods import WHAMAnalysis
        from pdynamo.us_data import USRun, read_pmf
        folder = self._us_folder(folder, sys, job)
        run = USRun(folder)
        files = sorted(glob.glob(os.path.join(folder, "data_collection", "window*.ptRes")))
        out_folder = os.path.join(folder, "wham")
        os.makedirs(out_folder, exist_ok=True)
        two_d = run.dims == 2
        parameters = {"analysis_type": "wham", "system": None, "file_names": files, "type": 1 if two_d else 0,
                      "logFrequency": 1, "maximumIterations": int(max_iterations),
                      "rmsGradientTolerance": float(rms_gradient),
                      "temperature": float(temperature or run.temperature),
                      "bins": [int(bins_rc2 or bins), int(bins)] if two_d else [int(bins)],
                      "resamples": int(resamples), "folder": out_folder, "logfile": name}
        ok, results = WHAMAnalysis().run(parameters)
        if not ok:
            raise APIError(str(results))
        dims, axes, F = read_pmf(results["pmf"])
        out = {"pmf_file": results["pmf"], "unit": "kJ/mol", "windows": len(files),
               "converged": "did not converge" not in results.get("msg", ""),
               "errors_file": results.get("pmf_errors")}
        if dims == 1:
            x = np.asarray(axes)
            out.update({"RC1": np.round(x, 4).tolist(), "PMF": np.round(F, 3).tolist()})
            out.update(_profile_1d(x, F))
            # . the histogram ends are thinly sampled: the PMF at the targets of
            #   the first and last windows is the better reference for the states
            targets = sorted(t[0] for t in (run.target(k) for k in run.keys) if t)
            ok = np.isfinite(F)
            if targets and ok.sum() > 1:
                at = lambda t: round(float(np.interp(t, x[ok], F[ok])), 3)
                out["first_window"] = {"rc": round(targets[0], 4), "energy": at(targets[0])}
                out["last_window"] = {"rc": round(targets[-1], 4), "energy": at(targets[-1])}
                if out.get("local_maxima"):
                    top = max(out["local_maxima"], key=lambda m: m["energy"])
                    out["barrier_from_first_window"] = round(top["energy"] - out["first_window"]["energy"], 3)
                    out["reaction_free_energy"] = round(out["last_window"]["energy"] - out["first_window"]["energy"], 3)
        else:
            xs, ys = axes
            out.update({"RC1": np.round(xs, 4).tolist(), "RC2": np.round(ys, 4).tolist(),
                        "PMF": [[None if not np.isfinite(v) else round(float(v), 3) for v in row] for row in F],
                        "note": "rows follow RC1, columns RC2"})
            out.update(_surface_2d(F, lambda i, j: xs[i], lambda i, j: ys[j], path_from, path_to))
        if self._summary(detail):
            _drop_arrays(out, ("RC1", "RC2", "PMF"), keep_short=False)
        return out

    def us_browser(self, folder=None, sys=None, job=None):
        """ Opens the Umbrella Sampling Browser on a run (to look at it with screenshot).

        Args:
          folder: the run folder (or sys and job of the umbrella sampling)
        """
        folder = self._us_folder(folder, sys, job)
        self.main.open_umbrella_browser(folder=folder)
        return {"opened": folder, "window": "Umbrella Sampling Browser"}

    # ============================================================ remote machines
    def remotes(self):
        """ Remote machines (Simulate > Remote Machines) and where simulations run now. """
        from remote import profiles as remote_profiles
        out = []
        for p in remote_profiles.load_all().values():
            out.append({"name": p["name"], "host": p["host"], "scheduler": p["scheduler"], "workdir": p["workdir"],
                        "cpus": p["cpus"], "partition": (p.get("slurm") or {}).get("partition") or None})
        return {"remotes": out, "run_on": getattr(self.p_session, "remote_target", None) or "this computer"}

    def remote_save(self, name, host, workdir, python="python3", setup="", scheduler="direct", cpus=1,
                    partition="", account="", time="24:00:00", mem="", fetch_exclude=None,
                    orca="", xtb="", dftb="", pdynamo_home="", sync_pdynamo=False):
        """ Adds or changes a remote machine (as Simulate > Remote Machines... does).

        Args:
          name: short name
          host: ssh destination: an alias of ~/.ssh/config, or user@host (keys only, no passwords)
          workdir: remote folder for the code and the jobs (on clusters, the scratch)
          python: interpreter that imports pDynamo after the setup lines
          setup: shell lines run first ("module load ...\nsource .../environment_bash.com")
          scheduler: "direct" (background process) or "slurm"
          partition, account, time, mem: SLURM options
          fetch_exclude: rsync patterns not fetched back, e.g. ["*.dcd"]
          orca, xtb, dftb: full paths of the QC programs there (parallel ORCA: load its MPI module in setup)
          pdynamo_home: the pDynamo folder there
          sync_pdynamo: jobs use a copy of it with this computer's pDynamo .py files and parameters
        """
        from remote import profiles as remote_profiles
        try:
            saved = remote_profiles.save({"name": name, "host": host, "workdir": workdir, "python": python,
                                          "setup": setup, "scheduler": scheduler, "cpus": int(cpus),
                                          "slurm": {"partition": partition, "account": account, "time": time, "mem": mem},
                                          "fetch_exclude": list(fetch_exclude or []),
                                          "qc_programs": {"orca": orca, "xtb": xtb, "dftb": dftb},
                                          "pdynamo_home": pdynamo_home, "sync_pdynamo": bool(sync_pdynamo)})
        except ValueError as error:
            raise APIError(str(error))
        try:
            self.main.refresh_run_on_menu()
        except Exception:
            pass
        return {"saved": saved["name"], "hint": "remote_test(name={!r})".format(saved["name"])}

    @off_main_thread
    def remote_test(self, name):
        """ Tests a remote machine: ssh, folder, scheduler, Python + pDynamo, and a real energy there.

        Args:
          name: the machine
        """
        from remote.jobs import RemoteJobs, RemoteError
        try:
            results = RemoteJobs().test(RemoteJobs().profile(name))
        except RemoteError as error:
            raise APIError(str(error))
        return {"ready": bool(results) and all(r["ok"] for r in results), "checks": results}

    @off_main_thread
    def remote_detect(self, host, name=None, save=False):
        """ Asks a machine (ssh) for what a remote profile needs and proposes one.

        Scheduler, partitions, accounts, scratch, the Python that really
        imports pDynamo, ORCA/xTB/DFTB+ (run to check them), the MPI module
        for parallel ORCA. Review the proposal; save=true stores it.

        Args:
          host: ssh alias of ~/.ssh/config, or user@host (the key must already work)
          name: profile name (default: the host)
          save: store the proposed profile (then remote_test)
        """
        from remote.discover import discover
        from remote import profiles as remote_profiles
        from remote import transport as remote_transport
        try:
            result = discover(host, name=name)
        except remote_transport.TransportError as error:
            raise APIError(str(error))
        if save:
            remote_profiles.save(result["profile"])
            result["saved"] = result["profile"]["name"]
        return result

    @off_main_thread
    def remote_partitions(self, name):
        """ The SLURM partitions of a remote machine (sinfo): time limit, nodes, idle CPUs, GPUs.

        Args:
          name: the machine
        """
        from remote import schedulers
        from remote import transport as remote_transport
        from remote.jobs import RemoteJobs, RemoteError
        try:
            profile = RemoteJobs().profile(name)
            return {"partitions": schedulers.get("slurm").partitions(remote_transport.make(profile))}
        except (RemoteError, RuntimeError, remote_transport.TransportError) as error:
            raise APIError(str(error))

    def run_on(self, remote=None):
        """ Where the next simulations run (Simulate > Run Simulations On).

        Args:
          remote: a machine name, or none / "local" for this computer
        """
        from remote import profiles as remote_profiles
        if remote in (None, "", "local", "this computer"):
            remote = None
        elif remote not in remote_profiles.load_all():
            raise APIError("No remote machine '{}' (see remotes).".format(remote))
        self.p_session.remote_target = remote
        try:
            self.main.refresh_run_on_menu()
        except Exception:
            pass
        return {"run_on": remote or "this computer"}

    @off_main_thread
    def remote_jobs(self, refresh=True, all=False):
        """ Remote jobs with their state (asks the machines first).

        Finished jobs are fetched automatically while EasyHybrid runs; use
        remote_fetch for a job of an earlier session.

        Args:
          refresh: ask the machines now (a few seconds)
          all: also the jobs that ended and were fetched
        """
        from remote.jobs import RemoteJobs, FINAL_STATES
        manager = RemoteJobs()
        records = manager.refresh() if refresh else list(manager.jobs().values())
        pending = [r["id"] for r in records if r.get("state") in ("finished", "error") and not r.get("fetched")]
        rows = getattr(self.p_session, "_remote_rows", {}) if self.p_session is not None else {}
        if pending and rows:
            self.p_session.remote_refresh_now([i for i in pending if i in rows])
        out = []
        for r in records:
            if not all and r.get("fetched") and r.get("state") in FINAL_STATES:
                continue
            out.append({k: r.get(k) for k in ("id", "name", "remote", "simulation_type", "state", "scheduler",
                                              "scheduler_id", "scheduler_state", "submitted", "local_folder",
                                              "logfile", "error", "fetched")})
            out[-1]["system_job"] = [r.get("system", {}).get("e_id"), r.get("system", {}).get("step")]
        return {"jobs": out}

    @off_main_thread
    def remote_fetch(self, remote_job):
        """ Copies the results of a remote job into its output folder now.

        Args:
          remote_job: remote job id (see remote_jobs)
        """
        from remote.jobs import RemoteJobs, RemoteError
        job = str(remote_job)
        try:
            record = RemoteJobs().fetch(job)
        except RemoteError as error:
            raise APIError(str(error))
        rows = getattr(self.p_session, "_remote_rows", {}) if self.p_session is not None else {}
        if job in rows:
            from gi.repository import GLib
            GLib.idle_add(self.p_session._remote_done, record)
        return {k: record.get(k) for k in ("id", "state", "local_folder", "logfile", "final_coordinates", "error")}

    def remote_cancel(self, remote_job):
        """ Cancels a remote job (scancel / kill).

        Args:
          remote_job: remote job id (see remote_jobs)
        """
        self.p_session.remote_cancel(str(remote_job))
        return {"cancelling": str(remote_job)}

    # ==================================================================== view
    def terminal(self, command):
        """ Runs an EasyHybrid terminal command (show, hide, color, center, select, zoom, dyff, opls, solvate...).

        "help" lists them, "help show" explains one. The reply is the text
        the terminal would print.

        Args:
          command: e.g. "show sticks obj=1 resn=HIS", "color obj=1 color=green symbol=C", "center obj=1 resi=10"
        """
        return {"output": self._cmd.run(str(command))}

    def set_visible(self, obj, visible=True, only=False):
        """ Shows/hides an object (the V checkbox of the Objects panel).

        Args:
          obj: object id
          visible: true/false
          only: hide every other object
        """
        target = self._vobject(obj)
        changed = []

        def walk(rows):
            for row in rows:
                yield row
                yield from walk(row.iterchildren())

        for row in walk(self.main.main_treeview.treestore):
            try:
                oid = int(row[1])
            except (TypeError, ValueError):
                continue
            if oid not in self.vm_session.vm_objects_dic or row.parent is None:
                continue
            want = bool(visible) if oid == target.index else (False if only else bool(row[6]))
            if bool(row[6]) != want:
                self.main.main_treeview.on_cell_visible_toggled(None, str(row.path))
                changed.append(oid)
        return {"changed": changed, "visible": [v.index for v in self.vm_session.vm_objects_dic.values() if v.active]}

    def set_frame(self, frame):
        """ Shows a frame (all objects; trajectories and scans).

        Args:
          frame: frame number (0-based)
        """
        self.vm_session.set_frame(int(frame))
        try:
            self.main.main_treeview.refresh_trajectory_scalebar()
        except Exception:
            pass
        return {"frame": int(frame)}

    def screenshot(self, path=None, target="viewer", scale=1, max_width=None):
        """ Saves a PNG of the 3D view or of a window (the AI can then look at it).

        Args:
          path: PNG file (default ~/.easyhybrid/shots/<time>.png)
          target: "viewer" (3D view only), "main" (main window) or part of a window title, e.g. "PES Analysis"
          scale: resolution multiple for the viewer (1-4)
          max_width: scale the PNG down to this width, 0 = full size (default: the AI settings)
        """
        if not path:
            shots = os.path.join(os.path.expanduser("~"), ".easyhybrid", "shots")
            os.makedirs(shots, exist_ok=True)
            path = os.path.join(shots, "{}_{}.png".format(time.strftime("%Y%m%d_%H%M%S"), target.replace(" ", "_")[:20]))
        path = os.path.abspath(os.path.expanduser(path))
        if target == "viewer":
            if not self.vm_session.vm_widget.save_image(path, scale_factor=int(scale)):
                raise APIError("Could not capture the 3D view.")
        else:
            from gi.repository import Gtk
            if target == "main":
                window = self.main.window
            else:
                matches = [w for w in Gtk.Window.list_toplevels()
                           if w.get_visible() and target.lower() in (w.get_title() or "").lower()]
                if not matches:
                    titles = [w.get_title() for w in Gtk.Window.list_toplevels() if w.get_visible() and w.get_title()]
                    raise APIError("No open window titled '{}'. Open: {}".format(target, titles))
                window = matches[0]
            self._draw_window_png(window, path)
        width = self._settings().get("screenshot_width", 0) if max_width is None else max_width
        size = None
        try:
            from PIL import Image
            image = Image.open(path)
            if width and image.width > int(width):
                height = round(image.height * int(width) / image.width)
                image = image.convert("RGB").resize((int(width), height), Image.LANCZOS)
                image.save(path)
            size = list(image.size)
        except ImportError:
            pass
        return {"path": path, "size": size}

    def _draw_window_png(self, window, path):
        """ PNG of a window drawn by GTK itself (no X11 grab: also Wayland/macOS);
            the OpenGL view, which GTK cannot draw, is pasted from its own capture. """
        import cairo
        import numpy as np
        alloc = window.get_allocation()
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, alloc.width, alloc.height)
        cr = cairo.Context(surface)
        window.draw(cr)
        gl = getattr(self.vm_session, "vm_widget", None)
        if gl is not None and gl.get_mapped() and gl.is_ancestor(window):
            image = gl.capture_screenshot()
            if image is not None:
                h, w = image.shape[:2]
                bgra = np.ascontiguousarray(image[..., [2, 1, 0, 3]])
                bgra[..., 3] = 255
                stride = cairo.ImageSurface.format_stride_for_width(cairo.FORMAT_ARGB32, w)
                if stride != 4 * w:
                    padded = np.zeros((h, stride // 4, 4), dtype=np.uint8)
                    padded[:, :w] = bgra
                    bgra = padded
                gl_surface = cairo.ImageSurface.create_for_data(memoryview(bgra), cairo.FORMAT_ARGB32, w, h, stride)
                x, y = gl.translate_coordinates(window, 0, 0)
                ga = gl.get_allocation()
                cr.save()
                cr.translate(x, y)
                cr.scale(ga.width / float(w), ga.height / float(h))
                cr.set_source_surface(gl_surface, 0, 0)
                cr.paint()
                cr.restore()
        surface.write_to_png(path)

    def windows(self):
        """ Titles of the open windows (for screenshot target= and dialogs). """
        from gi.repository import Gtk
        return {"windows": [w.get_title() for w in Gtk.Window.list_toplevels() if w.get_visible() and w.get_title()]}

    # ================================================================= dialogs
    def dialogs(self):
        """ Open dialogs (modal questions waiting for an answer), with their text and buttons. """
        try:
            from gi.repository import Gtk
        except ImportError:
            return []
        out = []
        for w in Gtk.Window.list_toplevels():
            if not isinstance(w, Gtk.Dialog) or not w.get_visible():
                continue
            texts, entries = [], []

            def walk(widget):
                if isinstance(widget, Gtk.Button):
                    return                      # . button labels are listed in "buttons"
                if isinstance(widget, Gtk.Label) and widget.get_visible():
                    t = widget.get_text().strip()
                    if t:
                        texts.append(t)
                if isinstance(widget, Gtk.Entry):
                    entries.append(widget.get_text())
                if isinstance(widget, Gtk.Container):
                    for child in widget.get_children():
                        walk(child)

            walk(w.get_content_area())
            buttons = []
            for b in w.get_action_area().get_children() if w.get_action_area() else []:
                if isinstance(b, Gtk.Button):
                    buttons.append({"label": (b.get_label() or "").replace("_", ""),
                                    "response": int(w.get_response_for_widget(b))})
            out.append({"id": id(w) & 0xFFFFFFFF, "title": w.get_title(), "modal": w.get_modal(),
                        "text": texts, "entries": entries, "buttons": buttons})
        return out

    def dialog_respond(self, button=None, dialog=None, response=None, text=None):
        """ Answers an open dialog (see dialogs): clicks a button by label or response id.

        Args:
          button: button label, e.g. "Yes", "No", "OK", "Import..."
          dialog: dialog id from `dialogs` (default: the most recent dialog)
          response: or the GTK response id
          text: text for the first entry of the dialog, if it has one
        """
        from gi.repository import Gtk
        dialogs = [w for w in Gtk.Window.list_toplevels() if isinstance(w, Gtk.Dialog) and w.get_visible()]
        if dialog is not None:
            dialogs = [w for w in dialogs if id(w) & 0xFFFFFFFF == int(dialog)]
        if not dialogs:
            raise APIError("No open dialog.")
        dialog = dialogs[-1]
        if text is not None:
            def first_entry(widget):
                if isinstance(widget, Gtk.Entry):
                    return widget
                if isinstance(widget, Gtk.Container):
                    for child in widget.get_children():
                        found = first_entry(child)
                        if found:
                            return found
                return None
            entry = first_entry(dialog.get_content_area())
            if entry is None:
                raise APIError("The dialog has no text entry.")
            entry.set_text(str(text))
        if response is None:
            if button is None:
                raise APIError("Give button= (label) or response=.")
            for b in dialog.get_action_area().get_children():
                if isinstance(b, Gtk.Button) and (b.get_label() or "").replace("_", "").lower() == str(button).lower():
                    response = dialog.get_response_for_widget(b)
                    break
            else:
                raise APIError("No button '{}' in the dialog.".format(button))
        dialog.response(int(response))
        return {"answered": dialog.get_title(), "response": int(response)}

    def events(self, since=0, kinds=None):
        """ Recorded events: finished jobs, errors, messages that dialogs would have shown.

        Args:
          since: only events with seq > since
          kinds: e.g. ["finished", "job_error", "error"]
        """
        from automation import events
        return {"events": events.since(int(since), kinds), "last": events.last_seq()}

    def python(self, code):
        """ Runs Python inside EasyHybrid (advanced; full access).

        Names available: main, p_session, vm_session, api, np. Assign the
        value to return to `result`; printed text is returned too.

        Args:
          code: Python source
        """
        import numpy as np
        namespace = {"main": self.main, "p_session": getattr(self.main, "p_session", None),
                     "vm_session": getattr(self.main, "vm_session", None), "api": self, "np": np, "result": None}
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            exec(compile(str(code), "<automation>", "exec"), namespace)
        return {"result": _jsonable(namespace.get("result")), "stdout": buffer.getvalue()[-20000:]}


# ====================================================================== schema
_PY_TYPES = {int: "integer", float: "number", bool: "boolean", str: "string", list: "array", dict: "object"}

# . arguments whose default is None: type and a default description
_KNOWN_ARGS = {
    "sys":       ({"type": ["integer", "string"]}, "system id (default: the active system)"),
    "obj":       ({"type": "integer"}, "object id (default: the last object of the system)"),
    "frame":     ({"type": "integer"}, "frame of the object (default: the last)"),
    "folder":    ({"type": "string"}, "output folder (default: the system's working folder)"),
    "atoms":     ({"type": ["array", "string"]}, "0-based atom indices: [0, 4, 7] or \"0-10,15\""),
    "of":        ({"type": ["array", "string"]}, "reference atoms for within="),
    "rc1":       ({"type": "object"}, None), "rc2": ({"type": "object"}, None),
    "parameters": ({"type": ["object", "array"]}, None),
    "job":       ({"type": "integer"}, "job number (default: the last job of the system)"),
    "name":      ({"type": "string"}, None), "path": ({"type": "string"}, None),
    "chain":     ({"type": "string"}, None), "resi": ({"type": ["string", "integer"]}, None),
    "resn":      ({"type": "string"}, None), "symbol": ({"type": "string"}, None),
    "within":    ({"type": "number"}, None), "free_within": ({"type": "number"}, None),
    "kinds":     ({"type": "array", "items": {"type": "string"}}, None),
    "reactant":  ({"type": "integer"}, None), "product": ({"type": "integer"}, None),
    "reactant_frame": ({"type": "integer"}, None), "product_frame": ({"type": "integer"}, None),
    "frames":    ({"type": "string"}, None), "logfile": ({"type": "string"}, None),
    "seed":      ({"type": "integer"}, None), "temperature_end": ({"type": "number"}, None),
    "pressure":  ({"type": "number"}, None), "response": ({"type": "integer"}, None),
    "dialog":    ({"type": "integer"}, None), "button": ({"type": "string"}, None),
    "text":      ({"type": "string"}, None), "topology": ({"type": "string"}, None),
    "parameter_folder": ({"type": "string"}, None), "scratch": ({"type": "string"}, None),
    "coordinates": ({"type": "string"}, None), "code": ({"type": "string"}, None),
    "command":   ({"type": "string"}, None),
    "wait":      ({"type": "boolean"}, "wait for the job and return its result in this call (default: the AI settings)"),
    "detail":    ({"type": "string", "enum": ["summary", "full"]}, "\"summary\" or \"full\" (default: the AI settings)"),
    "max_width": ({"type": "integer"}, None),
    "remote":    ({"type": "string"}, "remote machine name (see remotes); default: Simulate > Run Simulations On"),
    "remote_cpus": ({"type": "integer"}, None), "remote_job": ({"type": "string"}, None),
    "save": ({"type": "boolean"}, None),
    "host": ({"type": "string"}, None), "workdir": ({"type": "string"}, None), "fetch_exclude": ({"type": "array", "items": {"type": "string"}}, None),
    "source":    ({"type": "string"}, None), "force_constant": ({"type": "number"}, None),
    "bins_rc2":  ({"type": "integer"}, None), "temperature": ({"type": "number"}, None),
    "path_from": ({"type": "array", "items": {"type": "integer"}}, "2D only: start cell [i, j]"),
    "path_to":   ({"type": "array", "items": {"type": "integer"}}, "2D only: end cell [i, j]"),
}


def _arg_docs(doc):
    """ {"name": "description"} from the "Args:" block of a docstring;
        a line "a, b: text" documents both. """
    out, in_args, current = {}, False, None
    for line in (doc or "").splitlines():
        stripped = line.strip()
        if stripped.lower().startswith("args:"):
            in_args = True
            continue
        if not in_args:
            continue
        if not stripped:
            continue
        m = re.match(r"^([A-Za-z_][\w, *]*):\s*(.*)$", stripped)
        if m and len(line) - len(line.lstrip()) <= 10:
            names = [n.strip().rstrip("*").strip() for n in m.group(1).split(",")]
            current = names
            for n in names:
                out[n] = m.group(2)
        elif current:
            for n in current:
                out[n] += " " + stripped
    return out


def command_schemas():
    """ [{name, description, inputSchema}] of every public AutomationAPI method. """
    schemas = []
    for name, method in inspect.getmembers(AutomationAPI, predicate=inspect.isfunction):
        if name.startswith("_"):
            continue
        doc = inspect.getdoc(method) or ""
        summary = doc.split("\n\n")[0].replace("\n", " ").strip()
        body = doc.split("Args:")[0].strip()
        arg_docs = _arg_docs(doc)
        props, required = {}, []
        for pname, p in inspect.signature(method).parameters.items():
            if pname == "self":
                continue
            prop = {}
            default = p.default
            if default is inspect.Parameter.empty:
                required.append(pname)
            else:
                if type(default) in _PY_TYPES:
                    prop["type"] = _PY_TYPES[type(default)]
                if default is not None:
                    prop["default"] = default
            if "type" not in prop and pname in _KNOWN_ARGS:
                prop.update(_KNOWN_ARGS[pname][0])
            if pname in arg_docs:
                prop["description"] = arg_docs[pname]
            elif pname in _KNOWN_ARGS and _KNOWN_ARGS[pname][1]:
                prop["description"] = _KNOWN_ARGS[pname][1]
            props[pname] = prop
        schemas.append({"name": name, "summary": summary, "description": body,
                        "inputSchema": {"type": "object", "properties": props, "required": required}})
    return schemas
