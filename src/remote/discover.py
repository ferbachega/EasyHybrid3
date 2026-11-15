#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Detect Settings (2026-10-28): asks a machine reached by ssh what a remote
#  profile needs and proposes one -- scheduler, partitions and accounts,
#  a writable scratch folder, the Python + pDynamo that really import, the
#  QC programs (ORCA, xTB, DFTB+) and the MPI module for parallel ORCA.
#  Two round trips: an inventory script, then a probe of the combinations
#  found. Nothing is installed or changed on the machine.
#
#  discover(host) -> {"profile": {...}, "alternatives": {...}, "notes": [...],
#                     "inventory": {...}}

import os
import re
import subprocess

from remote import transport as remote_transport

_INVENTORY = r'''
echo "HOST	$(hostname)"
echo "USER	$USER"
echo "HOME	$HOME"
echo "NPROC	$(nproc 2>/dev/null)"
for c in sbatch squeue sacct sinfo scancel qsub bsub rsync python3; do
    p=$(command -v $c 2>/dev/null) && echo "CMD	$c	$p"
done
type module >/dev/null 2>&1 && echo "HAS_MODULE	1"
echo "SYSPY	$(python3 --version 2>&1)"
for d in "$SCRATCH" "/scratch/$USER" "/scratch/users/$USER" "/scratch1/$USER" "/work/$USER" "$WORK" \
         "/lustre/scratch/$USER" "/gpfs/scratch/$USER" "$HOME"; do
    [ -n "$d" ] && [ -d "$d" ] && [ -w "$d" ] && echo "SCRATCH	$d	$(df -Pk "$d" 2>/dev/null | awk 'NR==2{print $4}')"
done
if type module >/dev/null 2>&1; then
    module -t avail 2>&1 | grep -iE '^(python|anaconda|miniconda|conda|orca|xtb|dftb|openmpi)[^ ]*' | sed 's/^/MODULE	/'
    for m in $(module -t avail 2>&1 | grep -iE '^(orca|xtb|dftb)[^/]*/[^/ ]+' ); do
        for d in $(module show "$m" 2>&1 | grep -E 'prepend_path\("PATH"|prepend-path[[:space:]]+PATH' | grep -oE '/[^",)[:space:]]+'); do
            echo "MODPATH	$m	$d"
        done
    done
fi
for d in $(timeout 25 find "$HOME" -maxdepth 4 -type d -name pCore -path '*Dynamo*' 2>/dev/null | head -10); do
    d=${d%/pCore}; [ -d "$d/pBabel" ] && [ -d "$d/pMolecule" ] && [ -d "$d/parameters" ] && echo "PDYNAMO	$d"
done
for p in "$HOME"/programs/*orca*/orca "$HOME"/*orca*/orca "$HOME"/opt/*orca*/orca \
         "$HOME"/programs/xtb*/bin/xtb "$HOME"/xtb*/bin/xtb "$HOME"/opt/xtb*/bin/xtb \
         "$HOME"/programs/dftb*/bin/dftb+ "$HOME"/dftb*/bin/dftb+ "$HOME"/opt/dftb*/bin/dftb+; do
    [ -x "$p" ] && echo "PROG	$p"
done
for c in orca xtb dftb+; do p=$(command -v $c 2>/dev/null) && echo "PROG	$p"; done
for p in "$HOME"/miniconda3/envs/*/bin/dftb+ "$HOME"/anaconda3/envs/*/bin/dftb+ "$HOME"/miniforge3/envs/*/bin/dftb+ \
         "$HOME"/miniconda3/envs/*/bin/xtb "$HOME"/anaconda3/envs/*/bin/xtb "$HOME"/miniforge3/envs/*/bin/xtb; do
    [ -x "$p" ] && echo "PROG	$p"
done
for m in "$HOME"/opt/openmpi-*/bin/mpirun "$HOME"/programs/openmpi-*/bin/mpirun "$HOME"/openmpi-*/bin/mpirun; do
    [ -x "$m" ] && echo "MPIRUN	$m"
done
# . environment scripts the user keeps for these programs (what they export)
for f in "$HOME"/opt/*env*.sh "$HOME"/*env*.sh "$HOME"/programs/*env*.sh; do
    [ -f "$f" ] && echo "ENVFILE	$f	$(grep -oE 'export[[:space:]]+[A-Z0-9_]+' "$f" | awk '{print $2}' | sort -u | tr '\n' ' ')"
done
echo "WORKFREE	$(df -Pk "$HOME" 2>/dev/null | awk 'NR==2{print $4}')"
for e in "$HOME"/miniconda3/bin/python "$HOME"/anaconda3/bin/python "$HOME"/miniforge3/bin/python \
         "$HOME"/miniconda3/envs/*/bin/python "$HOME"/anaconda3/envs/*/bin/python "$HOME"/miniforge3/envs/*/bin/python; do
    [ -x "$e" ] && echo "CONDA	$e"
done
if command -v sinfo >/dev/null 2>&1; then
    out=$(sinfo -h -o '%P|%a|%l|%D|%C|%G' 2>&1); rc=$?
    if [ $rc -eq 0 ]; then echo "$out" | sed 's/^/PART	/'; else echo "SLURMERR	$(echo "$out" | head -1)"; fi
    sacctmgr -n -P show assoc user="$USER" format=account 2>/dev/null | sort -u | sed 's/^/ACCOUNT	/'
    sacctmgr -n -P show user "$USER" format=defaultaccount 2>/dev/null | sed 's/^/DEFACCOUNT	/'
fi
'''


def _parse(text):
    inv = {"cmd": {}, "scratch": [], "modules": [], "modpath": [], "pdynamo": [], "programs": [], "conda": [],
           "partitions": [], "accounts": [], "default_account": None}
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t")
        key, rest = parts[0], parts[1:]
        if not rest:
            continue
        if key in ("HOST", "USER", "HOME", "NPROC", "SYSPY"):
            inv[key.lower()] = rest[0]
        elif key == "CMD" and len(rest) == 2:
            inv["cmd"][rest[0]] = rest[1]
        elif key == "HAS_MODULE":
            inv["has_module"] = True
        elif key == "SCRATCH" and len(rest) == 2:
            inv["scratch"].append({"path": rest[0], "free_kb": int(rest[1]) if rest[1].isdigit() else None})
        elif key == "MODULE":
            name = rest[0].strip()
            if name and not name.endswith("/") and name not in inv["modules"]:
                inv["modules"].append(name)
        elif key == "MODPATH" and len(rest) == 2:
            inv["modpath"].append((rest[0], rest[1]))
        elif key == "PDYNAMO":
            inv["pdynamo"].append(rest[0])
        elif key == "PROG" and rest[0] not in inv["programs"]:
            inv["programs"].append(rest[0])
        elif key == "CONDA":
            inv["conda"].append(rest[0])
        elif key == "MPIRUN":
            inv.setdefault("mpirun", []).append(rest[0])
        elif key == "ENVFILE":
            inv.setdefault("envfiles", []).append((rest[0], (rest[1] if len(rest) > 1 else "").split()))
        elif key == "WORKFREE" and rest[0].isdigit():
            inv["home_free_kb"] = int(rest[0])
        elif key == "SLURMERR":
            inv["slurm_error"] = rest[0]
        elif key == "PART":
            from remote.schedulers import Slurm

            class _T:
                def run(self, *a, **k):
                    return 0, rest[0] + "\n", ""
            inv["partitions"] += Slurm().partitions(_T())
        elif key == "ACCOUNT" and rest[0].strip():
            inv["accounts"].append(rest[0].strip())
        elif key == "DEFACCOUNT" and rest[0].strip():
            inv["default_account"] = rest[0].strip()
    inv["partitions"].sort(key=lambda r: (not r["default"], r["name"]))
    return inv


def _version_key(text):
    return tuple(int(v) for v in re.findall(r"\d+", text)[:3])


def _python_modules(modules):
    """ python/3.x modules, newest first (3.8 or later). """
    out = []
    for m in modules:
        if not m.lower().startswith("python/"):
            continue
        v = _version_key(m.split("/", 1)[1])
        if len(v) >= 2 and v[0] == 3 and v[1] >= 8:
            out.append((v, m))
    return [m for _v, m in sorted(out, reverse=True)]


def _program_candidates(inv, name):
    """ Executables of one QC program: modules (also when the module's PATH is
        wrong, e.g. <prefix>/bin missing and the binary in <prefix>), the
        home folder and the PATH. [(path, origin)] """
    binary = {"orca": "orca", "xtb": "xtb", "dftb": "dftb+"}[name]
    found = []
    for module, path in inv["modpath"]:
        if module.lower().startswith(name):
            found.append((path.rstrip("/") + "/" + binary, "module " + module))
            parent = os.path.dirname(path.rstrip("/"))
            if path.rstrip("/").endswith("/bin"):
                found.append((parent + "/" + binary, "module {} (binary outside its bin/)".format(module)))
    for p in inv["programs"]:
        if os.path.basename(p) == binary:
            found.append((p, "home folder" if p.startswith(inv.get("home", "~")) else "PATH"))
    seen, out = set(), []
    for p, origin in found:
        if p not in seen:
            seen.add(p)
            out.append((p, origin))
    return out


def _probe_script(pdynamo_dirs, python_choices, programs, env_files=()):
    """ Second round trip: which Python imports which pDynamo, which program
        files are executable, the commit of each pDynamo. """
    lines = []
    for pd in pdynamo_dirs:
        env = pd + "/installation/shellScripts/environment_bash.com"
        lines.append('echo "COMMIT\t{0}\t$(git -C {0} rev-parse HEAD 2>/dev/null)"'.format(pd))
        lines.append('[ -f {0} ] && echo "ENVSCRIPT\t{1}\t{0}"'.format(env, pd))
        for kind, value in python_choices:
            load = "module purge >/dev/null 2>&1; module load {} >/dev/null 2>&1;".format(value) if kind == "module" else ""
            py = value if kind == "python" else "python3"
            lines.append(
                '( {load} [ -f {env} ] && source {env} >/dev/null 2>&1; export PYTHONPATH={pd}:$PYTHONPATH; '
                'r=$(timeout 60 {py} -c "import sys, numpy, yaml, pBabel, pMolecule; print(sys.version.split()[0])" 2>&1 | tail -1); '
                'echo "PROBE\t{pd}\t{kind}\t{value}\t$r" )'.format(load=load, env=env, pd=pd, py=py, kind=kind, value=value))
    for path in programs:
        lines.append('[ -e {0} ] && echo "REAL\t{0}\t$(readlink -f {0})"'.format(path))
    for path in programs:
        # . run it: a module binary may need the libraries its module loads
        lines.append('if [ -x {0} ] && [ -f {0} ]; then o=$(timeout 30 {0} --version 2>&1 | head -40); '
                     'if echo "$o" | grep -qiE "error while loading shared libraries|cannot execute|command not found"; '
                     'then echo "NORUN\t{0}\t$(echo "$o" | grep -iE "shared libraries|cannot execute|command not found" | head -1)"; '
                     'else echo "EXEC\t{0}\t$(echo "$o" | grep -oiE "(version|program version)[^0-9]*[0-9]+[.][0-9]+[.]?[0-9]*" | head -1)"; fi; fi'.format(path))
    for f in env_files:
        # . what mpirun a user script makes available (parallel ORCA)
        lines.append('( source {0} >/dev/null 2>&1; m=$(command -v mpirun) && echo "ENVMPI\t{0}\t$m\t$(mpirun --version 2>&1 | head -1)" )'.format(f))
    return "\n".join(lines) + "\n"


def _local_pdynamo_commit():
    home = os.environ.get("PDYNAMO3_HOME")
    if not home:
        return None
    try:
        return subprocess.run(["git", "-C", home, "rev-parse", "HEAD"], capture_output=True, text=True,
                              timeout=10).stdout.strip() or None
    except (OSError, subprocess.TimeoutExpired):
        return None


def discover(host, name=None, transport=None):
    """ Proposes a remote profile for host. Raises TransportError when ssh fails. """
    transport = transport or remote_transport.SSHTransport(host)
    notes = []
    rc, out, err = transport.run("bash -l -s", stdin=_INVENTORY, timeout=120)
    inv = _parse(out)
    if not inv.get("home"):
        raise remote_transport.TransportError("{}: no answer from the inventory ({}).".format(host, err.strip()[-300:]))

    # . scheduler
    if "sbatch" in inv["cmd"] and not inv.get("slurm_error"):
        scheduler = "slurm"
    elif "sbatch" in inv["cmd"]:
        scheduler = "direct"
        notes.append("SLURM is installed but does not answer ({}): jobs will run directly on {} "
                     "(background processes). Switch the scheduler to SLURM when its controller is up.".format(
                         inv["slurm_error"][:90], inv.get("host")))
    else:
        scheduler = "direct"
        if "qsub" in inv["cmd"] or "bsub" in inv["cmd"]:
            notes.append("This machine has {} (PBS/LSF), not supported yet: jobs will run as background processes "
                         "on the login node unless SLURM is available.".format("qsub" if "qsub" in inv["cmd"] else "bsub"))

    # . work folder: the first writable scratch that is not the home
    scratch = [s for s in inv["scratch"] if s["path"] != inv["home"]]
    base = scratch[0]["path"] if scratch else inv["home"]
    if not scratch:
        notes.append("No scratch file system found: the work folder is in the home folder.")

    # . Python x pDynamo
    pythons = [("system", "python3")] + [("module", m) for m in _python_modules(inv["modules"])[:3]] + \
              [("python", c) for c in inv["conda"][:3]]
    programs = {key: _program_candidates(inv, key) for key in ("orca", "xtb", "dftb")}
    env_files = [f for f, names in inv.get("envfiles", [])
                 if any(n in names for n in ("OMPI_HOME", "PDYNAMO3_ORCACOMMAND", "ORCA_DIR", "MPI_HOME"))]
    script = _probe_script(inv["pdynamo"][:3], pythons, [p for c in programs.values() for p, _o in c], env_files)
    rc, out, err = transport.run("bash -l -s", stdin=script, timeout=600)
    probes, executable, commits, envs, versions, norun, real, env_mpi = [], set(), {}, {}, {}, {}, {}, []
    for line in out.splitlines():
        parts = line.split("\t")
        if parts[0] == "PROBE" and len(parts) >= 5:
            ok = bool(re.match(r"^\d+\.\d+", parts[4].strip()))
            probes.append({"pdynamo": parts[1], "kind": parts[2], "value": parts[3], "ok": ok,
                           "detail": parts[4].strip()})
        elif parts[0] == "EXEC" and len(parts) >= 2:
            executable.add(parts[1])
            versions[parts[1]] = parts[2].strip() if len(parts) > 2 else ""
        elif parts[0] == "NORUN" and len(parts) >= 2:
            norun[parts[1]] = parts[2].strip() if len(parts) > 2 else "does not run"
        elif parts[0] == "COMMIT" and len(parts) == 3:
            commits[parts[1]] = parts[2].strip() or None
        elif parts[0] == "ENVSCRIPT" and len(parts) == 3:
            envs[parts[1]] = parts[2]
        elif parts[0] == "REAL" and len(parts) == 3:
            real[parts[1]] = parts[2].strip()
        elif parts[0] == "ENVMPI" and len(parts) >= 3:
            env_mpi.append((parts[1], parts[2], parts[3].strip() if len(parts) > 3 else ""))
    working = [p for p in probes if p["ok"]]

    setup, python, pdynamo_home = [], "python3", ""
    if working:
        # . prefer a module Python (newest), then the system one, then conda
        rank = {"module": 0, "system": 1, "python": 2}
        best = sorted(working, key=lambda p: (rank[p["kind"]], tuple(-v for v in _version_key(p["detail"]))))[0]
        pdynamo_home = best["pdynamo"]
        if best["kind"] == "module":
            setup.append("module load " + best["value"])
        elif best["kind"] == "python":
            python = best["value"]
        if pdynamo_home in envs:
            setup.append("source " + envs[pdynamo_home].replace(inv["home"], "~", 1))
        else:
            setup.append("export PYTHONPATH={}:$PYTHONPATH".format(pdynamo_home))
        failed = [p for p in probes if not p["ok"] and p["pdynamo"] == pdynamo_home]
        if failed:
            notes.append("Python {} imports pDynamo ({} {}); {} did not ({}).".format(
                best["detail"], best["kind"], best["value"],
                ", ".join("{} {}".format(p["kind"], p["value"]) for p in failed[:3]),
                "; ".join(p["detail"][:60] for p in failed[:2])))
    elif inv["pdynamo"]:
        notes.append("pDynamo was found ({}) but no Python imported it: {}.".format(
            ", ".join(inv["pdynamo"]), "; ".join("{} {}: {}".format(p["kind"], p["value"], p["detail"][:60])
                                                 for p in probes[:3])))
    else:
        notes.append("No pDynamo installation found in the home folder: install it there, or give its folder.")

    # . pDynamo of this computer: same commit -> send its Python files and parameters
    sync = False
    if pdynamo_home:
        local = _local_pdynamo_commit()
        remote_commit = commits.get(pdynamo_home)
        if local and remote_commit and local == remote_commit:
            sync = True
            notes.append("pDynamo there is at the same commit as here ({}): this computer's Python files and "
                         "parameters will be used there (local fixes included).".format(local[:7]))
        elif local and remote_commit:
            notes.append("pDynamo there is at commit {} and here at {}: not syncing this computer's files "
                         "(they could mismatch its compiled modules).".format(remote_commit[:7], local[:7]))

    # . QC programs
    qc_programs = {}
    for key, label in (("orca", "ORCA"), ("xtb", "xTB"), ("dftb", "DFTB+")):
        # . programs that run, the newest version first (from the version line, else from the path)
        # . one entry per real file (links), versioned before unversioned, newest first
        def version_of(po):
            return _version_key(versions.get(po[0]) or "") or _version_key(real.get(po[0], po[0])) or \
                _version_key(po[0])
        unique, seen_real = [], set()
        for p, o in programs[key]:
            if p in executable and real.get(p, p) not in seen_real:
                seen_real.add(real.get(p, p))
                unique.append((p, o))
        ok = sorted(unique, key=lambda po: (not version_of(po), [-v for v in version_of(po)]))
        bad_modules = [o for p, o in programs[key] if p not in executable and o.startswith("module") and "outside" not in o]
        if ok:
            qc_programs[key] = ok[0][0]
            if "outside" in ok[0][1]:
                notes.append("{}: {} adds a bin/ folder that does not hold the program; using {}.".format(
                    label, ok[0][1].split(" (")[0], ok[0][0]))
        elif bad_modules:
            notes.append("{}: {} found but no executable in its path.".format(label, ", ".join(bad_modules)))
        broken = [p for p, _o in programs[key] if p in norun]
        if broken:
            notes.append("{}: {} did not run alone ({}): it needs its module loaded.".format(
                label, ", ".join(broken), norun[broken[0]][:80]))
    if qc_programs.get("orca"):
        # . a user script that sets up ORCA with its MPI wins; then an Open MPI 4.1 module
        script_41 = [(f, m, v) for f, m, v in env_mpi if re.search(r"\b4\.1\.", v)]
        mpi = sorted((m for m in inv["modules"] if re.match(r"(?i)openmpi[/-]4\.1", m)),
                     key=lambda m: ("noncuda" not in m, "gcc" not in m, [-v for v in _version_key(m)]))
        if script_41:
            f, m, v = script_41[0]
            setup.append("source " + f.replace(inv["home"], "~", 1))
            notes.append("ORCA: {} sets up {} ({}); it was added to the setup lines for parallel ORCA.".format(
                f.replace(inv["home"], "~", 1), m.replace(inv["home"], "~", 1), v))
        elif mpi:
            setup.append("module load " + mpi[0])
            notes.append("Parallel ORCA (PAL) needs the MPI it was built with; {} was added to the setup "
                         "lines (ORCA 6 uses Open MPI 4.1).".format(mpi[0]))
        elif inv.get("mpirun"):
            notes.append("Open MPI found in {} but no script or module loads it: for parallel ORCA add its bin/ to "
                         "PATH and lib/ to LD_LIBRARY_PATH in the setup lines.".format(", ".join(inv["mpirun"])))
        else:
            notes.append("No Open MPI 4.1 found: parallel ORCA (PAL) may not run.")

    free_kb = (scratch[0]["free_kb"] if scratch else inv.get("home_free_kb"))
    if free_kb is not None and free_kb < 50 * 1024 * 1024:
        notes.append("Only {:.0f} GB free in the work folder: large runs (trajectories, umbrella sampling) "
                     "may fill it; use 'Do not fetch' and clean old jobs.".format(free_kb / 1024 / 1024))
    if pdynamo_home in envs:
        notes.append("The paths of the QC programs above override those of {} (it is still sourced for "
                     "PYTHONPATH and the pDynamo variables).".format(envs[pdynamo_home].replace(inv["home"], "~", 1)))

    # . SLURM
    slurm = {"partition": "", "account": "", "time": "24:00:00", "mem": "", "extra": ""}
    if scheduler == "slurm":
        if inv["partitions"]:
            slurm["partition"] = inv["partitions"][0]["name"]
        accounts = list(dict.fromkeys(inv["accounts"]))
        slurm["account"] = inv["default_account"] or (accounts[0] if accounts else "")

    profile = {"name": name or host, "host": host, "workdir": base.rstrip("/") + "/easyhybrid",
               "python": python, "setup": "\n".join(setup), "scheduler": scheduler, "cpus": 1,
               "slurm": slurm, "fetch_exclude": [], "qc_programs": {k: qc_programs.get(k, "") for k in ("orca", "xtb", "dftb")},
               "pdynamo_home": pdynamo_home.replace(inv["home"], "~", 1) if pdynamo_home else "",
               "sync_pdynamo": sync}
    alternatives = {
        "work folders": [s["path"] + "/easyhybrid" for s in inv["scratch"]],
        "python + pDynamo": ["{} {} -> {}".format(p["kind"], p["value"], p["pdynamo"]) for p in working],
        "partitions": [p["name"] for p in inv["partitions"]],
        "accounts": list(dict.fromkeys(inv["accounts"])),
        "orca": [p for p, _o in programs["orca"] if p in executable],
        "xtb": [p for p, _o in programs["xtb"] if p in executable],
        "dftb": [p for p, _o in programs["dftb"] if p in executable],
        "mpi modules": [m for m in inv["modules"] if m.lower().startswith("openmpi")],
        "mpi from scripts": ["{} -> {}".format(f, v) for f, _m, v in env_mpi],
    }
    return {"profile": profile, "alternatives": alternatives, "notes": notes,
            "inventory": {"hostname": inv.get("host"), "user": inv.get("user"), "cores": inv.get("nproc"),
                          "system python": inv.get("syspy"), "scheduler": scheduler,
                          "partitions": inv["partitions"]}}


def ssh_config_hosts(path=None):
    """ Host aliases of ~/.ssh/config (no wildcards). """
    path = path or os.path.join(os.path.expanduser("~"), ".ssh", "config")
    hosts = []
    try:
        for line in open(path):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0].lower() == "host":
                hosts += [h for h in parts[1:] if not any(c in h for c in "*?!")]
    except OSError:
        pass
    return list(dict.fromkeys(hosts))
