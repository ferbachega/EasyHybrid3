#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  ssh and rsync of the system (they honour ~/.ssh/config: aliases, keys,
#  ProxyJump). One shared connection per host (ControlMaster) makes the
#  frequent status checks cheap. BatchMode: a missing key or an unknown
#  host fails at once instead of waiting for a password nobody will type.
#
#  LocalTransport runs the same commands on this computer (tests).

import os
import shlex
import subprocess

CONTROL_DIR = os.path.join(os.path.expanduser("~"), ".easyhybrid")


class TransportError(Exception):
    pass


class SSHTransport:
    def __init__(self, host, connect_timeout=15):
        self.host = host
        self.connect_timeout = connect_timeout
        os.makedirs(CONTROL_DIR, mode=0o700, exist_ok=True)
        self.ssh_options = [
            "-o", "BatchMode=yes",
            "-o", "ConnectTimeout={}".format(connect_timeout),
            "-o", "ServerAliveInterval=30",
            "-o", "ControlMaster=auto",
            "-o", "ControlPath={}".format(os.path.join(CONTROL_DIR, "ssh-%C")),
            "-o", "ControlPersist=600",
        ]

    def run(self, command, timeout=120, check=True, stdin=None):
        """ Runs a shell command on the host: (returncode, stdout, stderr). """
        args = ["ssh"] + self.ssh_options + [self.host, command]
        try:
            p = subprocess.run(args, input=stdin, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise TransportError("ssh {} did not answer in {} s.".format(self.host, timeout))
        if check and p.returncode == 255:
            raise TransportError(_ssh_hint(self.host, p.stderr))
        return p.returncode, p.stdout, p.stderr

    def _remote(self, path):
        return "{}:{}".format(self.host, path)

    def rsync(self, sources, destination, upload=True, excludes=(), delete=False, timeout=3600, includes=()):
        """ rsync -a; sources/destination are local paths when not upload. With
            includes, only those patterns (then everything else excluded). """
        ssh = "ssh " + " ".join(shlex.quote(o) for o in self.ssh_options)
        args = ["rsync", "-a", "--partial", "-e", ssh]
        if delete:
            args.append("--delete")
        for pattern in excludes:
            args += ["--exclude", pattern]
        for pattern in includes:
            args += ["--include", pattern]
        if includes:
            args += ["--exclude", "*", "--prune-empty-dirs"]
        if upload:
            args += list(sources) + [self._remote(destination)]
        else:
            args += [self._remote(s) for s in sources] + [destination]
        try:
            p = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise TransportError("rsync with {} did not finish in {} s.".format(self.host, timeout))
        if p.returncode != 0:
            raise TransportError("rsync failed ({}): {}".format(p.returncode, (p.stderr or p.stdout).strip()[-800:]))

    def close(self):
        subprocess.run(["ssh"] + self.ssh_options + ["-O", "exit", self.host], capture_output=True, timeout=10)


class LocalTransport(SSHTransport):
    """ Same interface, on this computer (no ssh): for tests. """

    def __init__(self, host="localhost", connect_timeout=15):
        self.host = host
        self.ssh_options = []

    def run(self, command, timeout=120, check=True, stdin=None):
        p = subprocess.run(["bash", "-c", command], input=stdin, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr

    def rsync(self, sources, destination, upload=True, excludes=(), delete=False, timeout=3600, includes=()):
        args = ["rsync", "-a"] + (["--delete"] if delete else [])
        for pattern in excludes:
            args += ["--exclude", pattern]
        for pattern in includes:
            args += ["--include", pattern]
        if includes:
            args += ["--exclude", "*", "--prune-empty-dirs"]
        args += [os.path.expanduser(s) if not upload else s for s in sources] + [os.path.expanduser(destination)]
        p = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
        if p.returncode != 0:
            raise TransportError("rsync failed: {}".format(p.stderr.strip()))

    def close(self):
        pass


def _ssh_hint(host, stderr):
    text = (stderr or "").strip()
    low = text.lower()
    if "permission denied" in low:
        hint = "the ssh key was not accepted (EasyHybrid uses keys only; add yours with ssh-copy-id)"
    elif "host key verification failed" in low:
        hint = "the host key is unknown: connect once with 'ssh {}' in a terminal to accept it".format(host)
    elif "could not resolve" in low or "name or service not known" in low:
        hint = "the host name is unknown (an alias of ~/.ssh/config, or user@host)"
    elif "timed out" in low or "no route" in low or "unreachable" in low:
        hint = "the machine cannot be reached (VPN off? machine down?)"
    else:
        hint = "ssh failed"
    return "{}: {}. ({})".format(host, hint, text[-300:])


def make(profile):
    if profile.get("host") == "localhost-test":           # . tests: this computer, no ssh
        return LocalTransport()
    return SSHTransport(profile["host"])
