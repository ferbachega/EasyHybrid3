#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Remote jobs: bundle, submission, status, fetch, cancel, connection test.
#  The registry (~/.easyhybrid/remote_jobs.json) keeps every job, so the
#  results can be fetched after EasyHybrid was closed. GTK-free; the calls
#  that use ssh block, so the interface runs them in a thread.
#
#  A job record:
#    id, name, remote (profile name), host, scheduler, scheduler_id,
#    remote_dir (absolute), local_folder (where the results land, as a local
#    run would write them), simulation_type, system {label, tag, e_id, step},
#    submitted, state (submitted|queued|running|finished|error|failed|
#    cancelled), scheduler_state, logfile (local, once fetched), fetched,
#    error, warning

import copy
import json
import os
import pickle
import random
import shlex
import shutil
import string
import threading
import time

from remote import profiles as remote_profiles
from remote import schedulers
from remote import transport as remote_transport

CONTROL_DIR = os.path.join(os.path.expanduser("~"), ".easyhybrid")
REGISTRY_FILE = os.path.join(CONTROL_DIR, "remote_jobs.json")
SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE_TREES = ("pdynamo", "util", "remote")
CODE_EXCLUDES = ("__pycache__", "*.pyc", "old", "*.png")
FINAL_STATES = ("finished", "error", "failed", "cancelled")

_lock = threading.RLock()


class RemoteError(Exception):
    pass


def _q(text):
    return shlex.quote(str(text))


def qc_requirements(system, parameters):
    """ What the external QC program of a system needs remotely:
        {engine, threads, cpus, skf} -- cpus = parallel processes of the
        simulation x threads of the program (xTB "parallel", ORCA PALn). """
    import re
    model = getattr(system, "qcModel", None)
    name = type(model).__name__ if model is not None else ""
    engine = {"QCModelORCA": "ORCA", "QCModelXTB": "xTB", "QCModelDFTB": "DFTB+"}.get(name)
    threads = 1
    if engine == "xTB":
        threads = int(getattr(model, "parallel", 1) or 1)
    elif engine == "ORCA":
        text = " ".join(str(k) for k in (getattr(model, "keywords", None) or []))
        m = re.search(r"\bPAL(\d+)\b", text, re.I) or re.search(r"nprocs\s+(\d+)", text, re.I)
        threads = int(m.group(1)) if m else 1
    processes = int(parameters.get("NmaxThreads") or 1)
    skf = getattr(model, "skfPath", None) if engine == "DFTB+" else None
    return {"engine": engine, "threads": threads, "processes": processes, "cpus": processes * threads,
            "skf": skf if skf and os.path.isdir(skf) else None}


def _new_id():
    return time.strftime("%Y%m%d-%H%M%S-") + "".join(random.choices(string.ascii_lowercase + string.digits, k=4))


class RemoteJobs:
    def __init__(self, registry=None, profiles_file=None):
        self.registry = registry or REGISTRY_FILE
        self.profiles_file = profiles_file
        self._transports = {}

    # ------------------------------------------------------------ registry
    def jobs(self):
        with _lock:
            try:
                with open(self.registry) as f:
                    return {j["id"]: j for j in json.load(f).get("jobs", [])}
            except (OSError, ValueError):
                return {}

    def _save(self, records):
        with _lock:
            os.makedirs(os.path.dirname(self.registry), mode=0o700, exist_ok=True)
            tmp = self.registry + ".tmp"
            with open(tmp, "w") as f:
                json.dump({"jobs": list(records.values())}, f, indent=1, default=str)
            os.replace(tmp, self.registry)

    def _update(self, job_id, **changes):
        with _lock:
            records = self.jobs()
            if job_id not in records:
                raise RemoteError("No remote job {}.".format(job_id))
            records[job_id].update(changes)
            self._save(records)
            return records[job_id]

    def get(self, job_id):
        job = self.jobs().get(job_id)
        if job is None:
            raise RemoteError("No remote job {}.".format(job_id))
        return job

    def forget(self, job_id):
        with _lock:
            records = self.jobs()
            records.pop(job_id, None)
            self._save(records)

    # ------------------------------------------------------------ helpers
    def profile(self, name):
        try:
            return remote_profiles.get(name, self.profiles_file)
        except KeyError as error:
            raise RemoteError(str(error))

    def transport(self, profile):
        key = profile["host"]
        if key not in self._transports:
            self._transports[key] = remote_transport.make(profile)
        return self._transports[key]

    def _remote_root(self, profile, transport):
        """ The absolute work folder (created if needed). """
        workdir = profile["workdir"]
        if workdir.startswith("~/"):
            workdir = workdir[2:]
        rc, out, err = transport.run("mkdir -p {0} && cd {0} && pwd".format(_q(workdir)))
        if rc != 0 or not out.strip():
            raise RemoteError("Cannot use the work folder {} on {}: {}".format(workdir, profile["host"], err.strip()))
        return out.strip().splitlines()[-1]

    def sync_code(self, profile, transport=None, root=None):
        """ Copies src/pdynamo, src/util and src/remote to <workdir>/easyhybrid/src. """
        transport = transport or self.transport(profile)
        root = root or self._remote_root(profile, transport)
        code = root + "/easyhybrid/src"
        transport.run("mkdir -p {}".format(_q(code)))
        transport.rsync([os.path.join(SRC_DIR, t) for t in CODE_TREES], code + "/", upload=True,
                        excludes=CODE_EXCLUDES, delete=True)
        return code

    def sync_pdynamo(self, profile, transport=None, root=None):
        """ <workdir>/pdynamo3: a copy of the remote pDynamo (compiled modules
            kept) with this computer's Python files and parameters on top.
            Returns its path, or None when the profile does not ask for it. """
        if not profile.get("sync_pdynamo"):
            return None
        local = os.environ.get("PDYNAMO3_HOME")
        if not local or not os.path.isdir(os.path.join(local, "pCore")):
            raise RemoteError("PDYNAMO3_HOME is not set here: cannot send this computer's pDynamo files.")
        transport = transport or self.transport(profile)
        root = root or self._remote_root(profile, transport)
        target = root + "/pdynamo3"
        home = profile["pdynamo_home"]
        home = home[2:] if home.startswith("~/") else home
        rc, out, err = transport.run(
            "test -d {t}/pCore || rsync -a --exclude examples --exclude scratch {h}/ {t}/".format(t=_q(target), h=_q(home)),
            timeout=1800)
        if rc != 0:
            raise RemoteError("Could not copy pDynamo from {} on {}: {}".format(home, profile["host"], (err or out).strip()))
        transport.rsync([local + "/"], target + "/", upload=True, timeout=1800,
                        excludes=["examples/", "scratch/", "__pycache__/", "installation/", ".git/"],
                        includes=["*/", "*.py", "parameters/***"])
        return target

    # ------------------------------------------------------------ test
    def test(self, profile):
        """ Checks ssh, the work folder, the scheduler, Python + pDynamo, and
            runs a real energy on a small system. [ {check, ok, detail} ] """
        results = []

        def check(name, ok, detail=""):
            results.append({"check": name, "ok": bool(ok), "detail": str(detail).strip()[-600:]})
            return ok

        problems = remote_profiles.validate(profile)
        if not check("profile", not problems, "; ".join(problems)):
            return results
        transport = self.transport(profile)
        try:
            rc, out, err = transport.run("hostname", timeout=40)
        except remote_transport.TransportError as error:
            check("ssh", False, error)
            return results
        check("ssh", rc == 0, out or err)
        try:
            root = self._remote_root(profile, transport)
            check("work folder", True, root)
        except RemoteError as error:
            check("work folder", False, error)
            return results
        if profile["scheduler"] == "slurm":
            rc, out, err = transport.run("command -v sbatch squeue sacct scancel")
            check("SLURM commands", rc == 0 and out.count("\n") >= 4, out or err or "sbatch/squeue/sacct/scancel not found")
        try:
            code = self.sync_code(profile, transport, root)
            check("EasyHybrid code copied", True, code)
        except remote_transport.TransportError as error:
            check("EasyHybrid code copied", False, error)
            return results
        try:
            pdynamo_dir = self.sync_pdynamo(profile, transport, root)
            if pdynamo_dir:
                check("this computer's pDynamo files synced", True, pdynamo_dir)
        except (RemoteError, remote_transport.TransportError) as error:
            check("this computer's pDynamo files synced", False, error)
            return results
        # . external QC programs given in the profile
        for key, path in sorted((profile.get("qc_programs") or {}).items()):
            if path:
                rc, out, err = transport.run("test -x {0} && echo {0}".format(_q(path)))
                check("{} program".format({"orca": "ORCA", "xtb": "xTB", "dftb": "DFTB+"}.get(key, key)),
                      rc == 0, path if rc == 0 else "{} is not an executable file there".format(path))
        # . pDynamo + the code, with the profile's environment, on the login node
        env_pd = ("export PDYNAMO3_HOME={0}\nexport PDYNAMO3_PARAMETERS={0}/parameters\n"
                  "export PYTHONPATH={0}${{PYTHONPATH:+:$PYTHONPATH}}\n".format(_q(pdynamo_dir))) if pdynamo_dir else ""
        env_pd = env_pd.replace("{", "{{").replace("}", "}}")         # . formatted again below
        probe = ("{setup}\n" + env_pd + "export PYTHONPATH={code}${{PYTHONPATH:+:$PYTHONPATH}}\n"
                 "{py} -c 'import sys, numpy, pBabel, pMolecule, pdynamo.p_methods; "
                 "print(sys.version.split()[0], numpy.__version__)'").format(
            setup=profile.get("setup") or "", code=_q(code), py=profile.get("python") or "python3")
        rc, out, err = transport.run("bash -lc {}".format(_q(probe)), timeout=180)
        if not check("Python + pDynamo + EasyHybrid", rc == 0, (out if rc == 0 else err or out).strip().splitlines()[-1:] or ""):
            return results
        # . a real job: energy of a small MM system, run directly (not queued)
        try:
            job = self._smoke_job(profile, transport, root, code, env_pd)
            check("test energy (pickled system read and computed)", job[0], job[1])
        except Exception as error:
            check("test energy (pickled system read and computed)", False, error)
        # . the QC programs of the profile, through pDynamo there (catches a pDynamo that cannot read
        #   a program's output, missing libraries, wrong paths)
        for key, label in (("orca", "ORCA"), ("xtb", "xTB")):
            if (profile.get("qc_programs") or {}).get(key):
                try:
                    ok, detail = self._smoke_qc(profile, transport, root, code, env_pd, key)
                    check("{} energy through pDynamo there".format(label), ok, detail)
                except Exception as error:
                    check("{} energy through pDynamo there".format(label), False, error)
        return results

    _WATER = "3\nwater\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.469\nH 0.000 -0.757 -0.469\n"

    def _smoke_qc(self, profile, transport, root, code, env_pd, key):
        """ Energy of a water molecule with ORCA (HF/STO-3G) or xTB (GFN2) there;
            compared with this computer when the program exists here too. """
        from pBabel import ImportSystem, ExportSystem
        from pMolecule.QCModel import QCModelORCA, QCModelXTB
        tmp = os.path.join(CONTROL_DIR, "smoke_qc_" + key)
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
        xyz = os.path.join(tmp, "water.xyz")
        with open(xyz, "w") as f:
            f.write(self._WATER)
        system = ImportSystem(xyz)
        scratch = os.path.join(tmp, "scratch")
        os.makedirs(scratch)
        model = (QCModelORCA.WithOptions(keywords=["HF STO-3G"], scratch=scratch) if key == "orca" else
                 QCModelXTB.WithOptions(gfn=2, scratch=scratch))
        system.DefineQCModel(model)
        local = None
        try:
            local = system.Energy(log=None)
        except Exception:
            pass                                        # . the program may not exist here
        ExportSystem(os.path.join(tmp, "system.pkl"), system)
        remote_dir = root + "/jobs/_connection_test_" + key
        transport.run("rm -rf {0} && mkdir -p {0}".format(_q(remote_dir)))
        transport.rsync([os.path.join(tmp, "system.pkl")], remote_dir + "/", upload=True)
        from remote.profiles import QC_PROGRAM_VARIABLES
        exports = "".join("export {}={}\n".format(QC_PROGRAM_VARIABLES[k], _q(v))
                          for k, v in (profile.get("qc_programs") or {}).items() if v)
        script = ("{setup}\n" + env_pd + exports + "export PYTHONPATH={code}${{PYTHONPATH:+:$PYTHONPATH}}\ncd {dir}\n"
                  "{py} -c 'import os; from pBabel import ImportSystem; from remote.runner import _prepare_qc; "
                  "s = ImportSystem(\"system.pkl\"); _prepare_qc(s, os.getcwd()); print(\"ENERGY\", s.Energy(log=None))'").format(
            setup=profile.get("setup") or "", code=_q(code), dir=_q(remote_dir), py=profile.get("python") or "python3")
        rc, out, err = transport.run("bash -lc {}".format(_q(script)), timeout=300)
        lines = [l for l in out.splitlines() if l.startswith("ENERGY")]
        if rc != 0 or not lines or lines[-1].split()[1] in ("None", "nan"):
            tail = (err or out).strip().splitlines()
            return False, " | ".join(tail[-4:])[-500:] or "no energy"
        there = float(lines[-1].split()[1])
        if local is None:
            return True, "{:.4f} kJ/mol (the program is not available here to compare)".format(there)
        same = abs(there - local) < 0.05
        return same, "there {:.4f}, here {:.4f} kJ/mol{}".format(there, local, "" if same else
                                                                  ": different program versions or pDynamo?")

    def _smoke_job(self, profile, transport, root, code, env_pd=""):
        from pBabel import ImportSystem, ExportSystem
        from pMolecule.MMModel import MMModelOPLS
        from pMolecule.NBModel import NBModelFull
        examples = os.path.join(os.path.dirname(SRC_DIR), "examples", "tutorials", "01_bala", "bALA_not_opt.mol")
        system = ImportSystem(examples)
        system.DefineMMModel(MMModelOPLS.WithParameterSet("protein"))
        system.DefineNBModel(NBModelFull.WithDefaults())
        local_energy = system.Energy(log=None)
        tmp = os.path.join(CONTROL_DIR, "smoke")
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
        ExportSystem(os.path.join(tmp, "system.pkl"), system)
        remote_dir = root + "/jobs/_connection_test"
        transport.run("rm -rf {0} && mkdir -p {0}".format(_q(remote_dir)))
        transport.rsync([os.path.join(tmp, "system.pkl")], remote_dir + "/", upload=True)
        script = ("{setup}\n" + env_pd + "export PYTHONPATH={code}${{PYTHONPATH:+:$PYTHONPATH}}\ncd {dir}\n"
                  "{py} -c 'from pBabel import ImportSystem; s = ImportSystem(\"system.pkl\"); "
                  "print(\"ENERGY\", s.Energy(log=None))'").format(
            setup=profile.get("setup") or "", code=_q(code), dir=_q(remote_dir), py=profile.get("python") or "python3")
        rc, out, err = transport.run("bash -lc {}".format(_q(script)), timeout=300)
        lines = [l for l in out.splitlines() if l.startswith("ENERGY")]
        if rc != 0 or not lines:
            return False, (err or out).strip()[-600:]
        remote_energy = float(lines[-1].split()[1])
        same = abs(remote_energy - local_energy) < 1e-4 * max(1.0, abs(local_energy))
        return same, "remote {:.4f} kJ/mol, local {:.4f} kJ/mol{}".format(
            remote_energy, local_energy, "" if same else ": different pDynamo versions?")

    # ------------------------------------------------------------ submit
    def submit(self, system, parameters, profile_name, local_folder, cpus=None, name=None, overrides=None):
        """ prepare() + send(): blocks for the transfer. Returns the record. """
        job = self.prepare(system, parameters, profile_name, local_folder, cpus, name, overrides)
        return self.send(job["id"])

    def prepare(self, system, parameters, profile_name, local_folder, cpus=None, name=None, overrides=None):
        """ Writes the job bundle on this computer (no network; fast): the
            system, the parameters and the scheduler script. The record is
            saved with state "preparing"; send() transfers and submits it.
            parameters: as run_simulation() gets them, "system" left out. """
        profile = self.profile(profile_name)
        scheduler = schedulers.get(profile["scheduler"])
        job_id = _new_id()
        sim_type = parameters["simulation_type"]
        name = name or parameters.get("traj_folder_name") or parameters.get("trajectory_name") or \
            parameters.get("vobject_name") or sim_type
        qc = qc_requirements(system, parameters)
        bundle = os.path.join(CONTROL_DIR, "bundles", job_id)
        os.makedirs(bundle)
        try:
            self._write_system(system, os.path.join(bundle, "system.pkl"))
            if qc["skf"]:
                # . DFTB+ Slater-Koster files travel with the job (runner: input/skf)
                shutil.copytree(qc["skf"], os.path.join(bundle, "input", "skf"))
            params = {k: v for k, v in parameters.items() if k not in ("system", "queue", "logfile")}
            params["folder"] = "output"
            source = params.get("source_folder")
            if source:
                params["source_folder"] = "input/" + os.path.basename(os.path.normpath(source))
            with open(os.path.join(bundle, "parameters.pkl"), "wb") as f:
                pickle.dump(params, f, protocol=4)
        except Exception:
            shutil.rmtree(bundle, ignore_errors=True)
            raise
        job = {"id": job_id, "name": str(name), "remote": profile["name"], "host": profile["host"],
               "scheduler": scheduler.name, "local_folder": os.path.abspath(local_folder),
               "simulation_type": sim_type, "source_folder": os.path.normpath(source) if source else None,
               "cpus": int(cpus or (qc["cpus"] if qc["engine"] or qc["processes"] > 1 else 0) or profile.get("cpus") or 1),
               "qc_engine": qc["engine"], "qc_threads": qc["threads"], "processes": qc["processes"],
               "system": {"label": getattr(system, "label", None), "tag": getattr(system, "e_tag", None),
                          "e_id": getattr(system, "e_id", None), "step": getattr(system, "e_step_counter", None)},
               "submitted": time.strftime("%Y-%m-%d %H:%M:%S"), "state": "preparing", "fetched": False,
               "bundle": bundle}
        job.update(overrides or {})
        with _lock:
            records = self.jobs()
            records[job_id] = job
            self._save(records)
        return job

    def send(self, job_id):
        """ Transfers a prepared bundle (and the code) and submits it. On
            failure the record gets state "error" and the error is raised. """
        job = self.get(job_id)
        profile = self.profile(job["remote"])
        scheduler = schedulers.get(job["scheduler"])
        transport = self.transport(profile)
        bundle = job["bundle"]
        try:
            root = self._remote_root(profile, transport)
            code = self.sync_code(profile, transport, root)
            pdynamo_dir = self.sync_pdynamo(profile, transport, root)
            if pdynamo_dir:
                profile = dict(profile, _pdynamo_dir=pdynamo_dir)
            remote_dir = "{}/jobs/{}".format(root, job_id)
            job["remote_dir"] = remote_dir
            record = {k: v for k, v in job.items() if k != "bundle"}
            with open(os.path.join(bundle, "job.json"), "w") as f:
                json.dump(record, f, indent=1, default=str)
            with open(os.path.join(bundle, scheduler.script_name), "w") as f:
                f.write(scheduler.script(profile, job, remote_dir, code))
            transport.run("mkdir -p {}".format(_q(root + "/jobs")))
            transport.rsync([bundle], root + "/jobs/", upload=True)
            if job.get("source_folder"):
                transport.run("mkdir -p {}".format(_q(remote_dir + "/input")))
                transport.rsync([job["source_folder"]], remote_dir + "/input/", upload=True)
            scheduler_id = scheduler.submit(transport, remote_dir)
        except Exception as error:
            self._update(job_id, state="error", error="Not submitted: {}".format(error))
            raise
        finally:
            shutil.rmtree(bundle, ignore_errors=True)
        return self._update(job_id, remote_dir=remote_dir, scheduler_id=scheduler_id, bundle=None,
                            state="queued" if scheduler.name == "slurm" else "running")

    @staticmethod
    def _write_system(system, path):
        """ ExportSystem without the interface objects and the session history. """
        from pBabel import ExportSystem
        saved = {}
        for attr, empty in (("e_treeview_iter", None), ("e_liststore_iter", None),
                            ("e_job_history", {}), ("e_logfile_data", {})):
            if hasattr(system, attr):
                saved[attr] = getattr(system, attr)
                setattr(system, attr, copy.copy(empty))
        try:
            ExportSystem(path, system)
        finally:
            for attr, value in saved.items():
                setattr(system, attr, value)

    # ------------------------------------------------------------ status
    def refresh(self, job_ids=None):
        """ Asks the machines about the unfinished jobs; returns the records. """
        records = self.jobs()
        todo = [j for j in records.values() if j.get("state") not in FINAL_STATES
                and j.get("scheduler_id") and (job_ids is None or j["id"] in job_ids)]
        by_remote = {}
        for job in todo:
            by_remote.setdefault(job["remote"], []).append(job)
        for remote_name, jobs in by_remote.items():
            try:
                profile = self.profile(remote_name)
                transport = self.transport(profile)
                scheduler = schedulers.get(jobs[0]["scheduler"])
                states = scheduler.states(transport, [j["scheduler_id"] for j in jobs])
                statuses = self._read_status(transport, [j["remote_dir"] for j in jobs])
            except (RemoteError, remote_transport.TransportError, KeyError) as error:
                for job in jobs:
                    self._update(job["id"], last_error=str(error), checked=time.strftime("%H:%M:%S"))
                continue
            for job in jobs:
                sched = states.get(str(job["scheduler_id"]), "GONE")
                status = statuses.get(job["remote_dir"]) or {}
                if status.get("state") in ("finished", "error"):
                    state = status["state"]
                elif sched in schedulers.QUEUED:
                    state = "queued"
                elif sched in schedulers.ACTIVE or sched == "RUNNING":
                    state = "running"
                elif sched.startswith("CANCELLED"):
                    state = "cancelled"
                else:
                    state = "failed"          # . the process or the SLURM job ended without a status
                changes = {"state": state, "scheduler_state": sched, "checked": time.strftime("%H:%M:%S"),
                           "last_error": None}
                if status.get("error"):
                    changes["error"] = status["error"]
                if status.get("warning"):
                    changes["warning"] = status["warning"]
                if status.get("logfile"):
                    changes["remote_logfile"] = status["logfile"]
                self._update(job["id"], **changes)
        records = self.jobs()
        return [records[j] for j in (job_ids or records)] if job_ids else list(records.values())

    @staticmethod
    def _read_status(transport, dirs):
        marker = "@@EASYHYBRID@@"
        cmd = "; ".join("echo {m} {d}; cat {d}/status.json 2>/dev/null".format(m=marker, d=_q(d)) for d in dirs)
        _, out, _ = transport.run(cmd)
        result, current, buf = {}, None, []
        for line in out.splitlines() + [marker + " END"]:
            if line.startswith(marker):
                if current is not None:
                    try:
                        result[current] = json.loads("\n".join(buf)) if buf else {}
                    except ValueError:
                        result[current] = {}
                current, buf = line[len(marker):].strip(), []
            else:
                buf.append(line)
        return result

    # ------------------------------------------------------------ fetch / cancel
    def fetch(self, job_id, exclude=None):
        """ Copies output/ of the job into its local folder (as a local run
            would have written it) and the logs into
            <local folder>/.easyhybrid_remote/<id>/. Returns the record. """
        job = self.get(job_id)
        profile = self.profile(job["remote"])
        transport = self.transport(profile)
        os.makedirs(job["local_folder"], exist_ok=True)
        excludes = list(profile.get("fetch_exclude") or []) + list(exclude or [])
        transport.rsync([job["remote_dir"] + "/output/"], job["local_folder"] + "/", upload=False, excludes=excludes)
        logs = os.path.join(job["local_folder"], ".easyhybrid_remote", job_id)
        os.makedirs(logs, exist_ok=True)
        for name in ("status.json", "runner.log", "scheduler.log"):
            try:
                transport.rsync([job["remote_dir"] + "/" + name], logs + "/", upload=False)
            except remote_transport.TransportError:
                pass
        status = {}
        try:
            with open(os.path.join(logs, "status.json")) as f:
                status = json.load(f)
        except (OSError, ValueError):
            pass
        changes = {"fetched": True, "fetched_at": time.strftime("%Y-%m-%d %H:%M:%S"), "logs": logs}
        if status.get("state") in ("finished", "error"):
            changes["state"] = status["state"]
        if status.get("logfile"):
            changes["logfile"] = os.path.join(job["local_folder"], status["logfile"])
        if status.get("final_coordinates"):
            changes["final_coordinates"] = os.path.join(job["local_folder"], status["final_coordinates"])
        if status.get("error"):
            changes["error"] = status["error"]
        return self._update(job_id, **changes)

    def cancel(self, job_id):
        job = self.get(job_id)
        profile = self.profile(job["remote"])
        schedulers.get(job["scheduler"]).cancel(self.transport(profile), job["scheduler_id"])
        return self._update(job_id, state="cancelled")

    def clean(self, job_id):
        """ Deletes the job folder on the remote machine (after fetching). """
        job = self.get(job_id)
        if not job.get("remote_dir", "").rstrip("/").split("/")[-1] == job_id:
            raise RemoteError("Refusing to delete an unexpected folder: {}".format(job.get("remote_dir")))
        profile = self.profile(job["remote"])
        self.transport(profile).run("rm -rf {}".format(_q(job["remote_dir"])))
        return self._update(job_id, cleaned=True)
