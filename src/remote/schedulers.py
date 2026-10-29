#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  How a job is started, followed and cancelled on the remote machine.
#  The job folder holds the script written here; the runner (runner.py)
#  writes status.json, which says whether the simulation itself succeeded.
#  The scheduler only says queued / running / gone.
#
#    direct  a background process (nohup setsid), for workstations
#    slurm   sbatch / squeue / sacct / scancel
#
#  A new scheduler: a class with script(), submit(), states(), cancel()
#  and an entry in SCHEDULERS.

import re
import shlex


def _q(text):
    return shlex.quote(str(text))


def _body(profile, job_dir, code_dir, job=None):
    """ The shell lines shared by every scheduler: environment, then the runner. """
    job = job or {}
    lines = ["set -o pipefail", "cd {}".format(_q(job_dir))]
    setup = (profile.get("setup") or "").strip()
    if setup:
        lines += ["# . environment (profile 'setup')", setup]
    from remote.profiles import QC_PROGRAM_VARIABLES
    for key, variable in QC_PROGRAM_VARIABLES.items():
        path = ((profile.get("qc_programs") or {}).get(key) or "").strip()
        if path:
            lines.append("export {}={}".format(variable, _q(path)))
    pdynamo_dir = profile.get("_pdynamo_dir")
    if pdynamo_dir:
        # . the synced pDynamo copy wins over the PYTHONPATH of the setup lines
        lines += ["export PDYNAMO3_HOME={}".format(_q(pdynamo_dir)),
                  "export PDYNAMO3_PARAMETERS=$PDYNAMO3_HOME/parameters",
                  "export PDYNAMO3_STYLE=$PDYNAMO3_PARAMETERS/ccsStyleSheets/defaultStyle.css",
                  "export PYTHONPATH=$PDYNAMO3_HOME${PYTHONPATH:+:$PYTHONPATH}"]
    lines += ["export PYTHONPATH={}${{PYTHONPATH:+:$PYTHONPATH}}".format(_q(code_dir)),
              "export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}",
              # . several ORCA mpiruns at once (2D scan, umbrella windows): no core binding,
              #   or they all pile up on the same cores
              ("export OMPI_MCA_hwloc_base_binding_policy=none" if job.get("qc_engine") == "ORCA"
               and int(job.get("processes") or 1) > 1 else "# ."),
              "{} -u -m remote.runner {} > runner.log 2>&1".format(profile.get("python") or "python3", _q(job_dir))]
    return "\n".join(lines) + "\n"


class Direct:
    """ A background process on the remote machine (no queue). """
    name = "direct"
    script_name = "run.sh"

    def script(self, profile, job, job_dir, code_dir):
        return "#!/bin/bash -l\n# EasyHybrid job {}\n".format(job["id"]) + _body(profile, job_dir, code_dir, job)

    def submit(self, transport, job_dir):
        # . bash -l: the profile's 'setup' may use "module", defined by login shells
        cmd = "cd {0} && nohup setsid bash -l {1} > scheduler.log 2>&1 < /dev/null & echo $!".format(_q(job_dir), self.script_name)
        rc, out, err = transport.run(cmd)
        pid = out.strip().splitlines()[-1] if out.strip() else ""
        if rc != 0 or not pid.isdigit():
            raise RuntimeError("Could not start the job: {}".format((err or out).strip()))
        return pid

    def states(self, transport, ids):
        """ {id: "RUNNING" or "GONE"} """
        if not ids:
            return {}
        cmd = "for p in {}; do if kill -0 $p 2>/dev/null; then echo $p RUNNING; else echo $p GONE; fi; done".format(
            " ".join(_q(i) for i in ids))
        _, out, _ = transport.run(cmd)
        return dict(line.split() for line in out.splitlines() if len(line.split()) == 2)

    def cancel(self, transport, job_id):
        # . the whole process group (setsid made the job its leader)
        transport.run("kill -TERM -- -{0} 2>/dev/null || kill -TERM {0}".format(_q(job_id)), check=False)


class Slurm:
    name = "slurm"
    script_name = "job.slurm"

    def script(self, profile, job, job_dir, code_dir):
        s = dict(profile.get("slurm") or {})
        cpus = int(job.get("cpus") or profile.get("cpus") or 1)
        # . parallel ORCA runs mpirun, and Open MPI counts SLURM tasks (slots), not CPUs
        mpi = job.get("qc_engine") == "ORCA" and int(job.get("qc_threads") or 1) > 1
        lines = ["#!/bin/bash -l",
                 "#SBATCH --job-name=eh_{}".format(re.sub(r"[^A-Za-z0-9_.-]", "_", job.get("name") or job["id"])[:40]),
                 "#SBATCH --nodes=1",
                 "#SBATCH --ntasks={}".format(cpus if mpi else 1),
                 "#SBATCH --cpus-per-task={}".format(1 if mpi else cpus),
                 "#SBATCH --output=scheduler.log"]
        for key, flag in (("partition", "--partition"), ("account", "--account"), ("time", "--time"), ("mem", "--mem")):
            value = job.get(key) or s.get(key)
            if value:
                lines.append("#SBATCH {}={}".format(flag, value))
        for extra in (s.get("extra") or "").splitlines():
            if extra.strip():
                lines.append(extra.strip() if extra.strip().startswith("#SBATCH") else "#SBATCH " + extra.strip())
        lines.append("# EasyHybrid job {}".format(job["id"]))
        return "\n".join(lines) + "\n" + _body(profile, job_dir, code_dir, job)

    def submit(self, transport, job_dir):
        rc, out, err = transport.run("cd {} && sbatch --parsable {}".format(_q(job_dir), self.script_name))
        job_id = out.strip().split(";")[0] if out.strip() else ""
        if rc != 0 or not job_id.isdigit():
            raise RuntimeError("sbatch failed: {}".format((err or out).strip()))
        return job_id

    def states(self, transport, ids):
        """ {id: SLURM state}: PENDING, RUNNING, COMPLETED, FAILED, CANCELLED, TIMEOUT, ... """
        if not ids:
            return {}
        joined = ",".join(str(i) for i in ids)
        out_states = {}
        _, out, _ = transport.run("squeue -h -j {} -o '%i %T' 2>/dev/null".format(joined), check=True)
        for line in out.splitlines():
            parts = line.split()
            if len(parts) == 2:
                out_states[parts[0]] = parts[1]
        missing = [i for i in ids if str(i) not in out_states]
        if missing:
            _, out, _ = transport.run("sacct -n -X -P -j {} -o JobID,State 2>/dev/null".format(",".join(missing)))
            for line in out.splitlines():
                parts = line.split("|")
                if len(parts) == 2:
                    out_states[parts[0]] = parts[1].split()[0]
        for i in ids:
            out_states.setdefault(str(i), "GONE")
        return out_states

    def cancel(self, transport, job_id):
        transport.run("scancel {}".format(_q(job_id)), check=False)

    def partitions(self, transport):
        """ The partitions of the cluster (sinfo): [{name, default, available,
            time_limit, nodes, cpus_idle, cpus_total, gres}], default first. """
        rc, out, err = transport.run("sinfo -h -o '%P|%a|%l|%D|%C|%G'")
        if rc != 0:
            raise RuntimeError("sinfo failed: {}".format((err or out).strip()))
        rows = []
        for line in out.splitlines():
            parts = line.strip().split("|")
            if len(parts) < 6:
                continue
            name, avail, limit, nodes, cpus, gres = parts[:6]
            try:
                _alloc, idle, _other, total = (int(v) for v in cpus.split("/"))
            except ValueError:
                idle = total = None
            rows.append({"name": name.rstrip("*"), "default": name.endswith("*"), "available": avail,
                         "time_limit": limit, "nodes": int(nodes) if nodes.isdigit() else None,
                         "cpus_idle": idle, "cpus_total": total,
                         "gres": None if gres in ("(null)", "") else gres})
        rows.sort(key=lambda r: (not r["default"], r["name"]))
        return rows


SCHEDULERS = {"direct": Direct(), "slurm": Slurm()}

QUEUED = {"PENDING", "CONFIGURING", "REQUEUED", "SUSPENDED", "RESV_DEL_HOLD", "REQUEUE_HOLD"}
ACTIVE = {"RUNNING", "COMPLETING", "STAGE_OUT"}


def get(name):
    if name not in SCHEDULERS:
        raise KeyError("Unknown scheduler '{}' ({}).".format(name, ", ".join(SCHEDULERS)))
    return SCHEDULERS[name]
