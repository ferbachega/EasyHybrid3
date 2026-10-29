#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Remote machine profiles, saved in ~/.easyhybrid/remotes.json. No password
#  is stored: authentication is by ssh key, and "host" is best an alias of
#  ~/.ssh/config (user, key, port and ProxyJump come from there).

import json
import os

CONTROL_DIR = os.path.join(os.path.expanduser("~"), ".easyhybrid")
PROFILES_FILE = os.path.join(CONTROL_DIR, "remotes.json")

SCHEDULERS = ("direct", "slurm")

DEFAULT_PROFILE = {
    "name": "",
    "host": "",                 # ssh destination: alias of ~/.ssh/config, or user@host
    "workdir": "",              # remote folder for the code and the jobs (absolute, or relative to $HOME)
    "python": "python3",        # interpreter that imports pDynamo
    "setup": "",                # shell lines run before python: module load, source pDynamo env, conda...
    "scheduler": "direct",      # direct (background process) or slurm
    "cpus": 1,                  # default CPUs of a job (2D scans and umbrella sampling use them)
    "slurm": {"partition": "", "account": "", "time": "24:00:00", "mem": "", "extra": ""},
    "fetch_exclude": [],        # rsync patterns left on the remote machine when fetching, e.g. "*.dcd"
    # . external QC programs on the remote machine (empty: whatever the setup lines define);
    #   exported as PDYNAMO3_ORCACOMMAND / PDYNAMO3_XTBCOMMAND / PDYNAMO3_DFTBCOMMAND
    "qc_programs": {"orca": "", "xtb": "", "dftb": ""},
    # . pDynamo: the remote installation (its compiled modules suit the remote Python); with
    #   sync_pdynamo, jobs use a copy of it in <workdir>/pdynamo3 with this computer's Python
    #   files and parameters on top (local fixes and DYFF parameters), the original untouched
    "pdynamo_home": "",
    "sync_pdynamo": False,
}

QC_PROGRAM_VARIABLES = {"orca": "PDYNAMO3_ORCACOMMAND", "xtb": "PDYNAMO3_XTBCOMMAND", "dftb": "PDYNAMO3_DFTBCOMMAND"}


def _merge(profile):
    out = json.loads(json.dumps(DEFAULT_PROFILE))
    for key, value in (profile or {}).items():
        if key in ("slurm", "qc_programs") and isinstance(value, dict):
            out[key].update(value)
        elif key in out:
            out[key] = value
    return out


def load_all(path=None):
    """ {name: profile} """
    try:
        with open(path or PROFILES_FILE) as f:
            data = json.load(f)
    except (OSError, ValueError):
        return {}
    return {p["name"]: _merge(p) for p in data.get("profiles", []) if p.get("name")}


def get(name, path=None):
    profiles = load_all(path)
    if name not in profiles:
        raise KeyError("No remote machine named '{}' (known: {}).".format(name, ", ".join(profiles) or "none"))
    return profiles[name]


def validate(profile):
    """ List of problems (empty when the profile is usable). """
    problems = []
    if not profile.get("name"):
        problems.append("the profile has no name")
    if not profile.get("host"):
        problems.append("no ssh host")
    if not profile.get("workdir"):
        problems.append("no remote work folder")
    if profile.get("sync_pdynamo") and not profile.get("pdynamo_home"):
        problems.append("syncing pDynamo needs the folder of pDynamo on the remote machine")
    if profile.get("scheduler") not in SCHEDULERS:
        problems.append("scheduler must be one of {}".format(", ".join(SCHEDULERS)))
    try:
        if int(profile.get("cpus", 1)) < 1:
            problems.append("cpus must be at least 1")
    except (TypeError, ValueError):
        problems.append("cpus must be a number")
    return problems


def save(profile, path=None):
    """ Adds or replaces a profile (by name). """
    profile = _merge(profile)
    problems = validate(profile)
    if problems:
        raise ValueError("; ".join(problems))
    path = path or PROFILES_FILE
    profiles = load_all(path)
    profiles[profile["name"]] = profile
    os.makedirs(os.path.dirname(path), mode=0o700, exist_ok=True)
    with open(path, "w") as f:
        json.dump({"profiles": list(profiles.values())}, f, indent=1)
    return profile


def remove(name, path=None):
    path = path or PROFILES_FILE
    profiles = load_all(path)
    if profiles.pop(name, None) is None:
        return False
    with open(path, "w") as f:
        json.dump({"profiles": list(profiles.values())}, f, indent=1)
    return True
