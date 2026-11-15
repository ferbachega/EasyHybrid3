#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Remote Machines > Add SSH Host... (2026-10-28): writes a "Host" block to
#  ~/.ssh/config, so a machine given as address + user + port + key becomes
#  an alias that ssh, rsync and EasyHybrid all use. Never touches an
#  existing alias; keeps a backup; no passwords anywhere.
#
#  The two steps that need the password of the remote account (authorizing
#  the key) are given as commands for the user's own terminal:
#  setup_commands().

import os
import re
import shutil
import time

SSH_DIR = os.path.join(os.path.expanduser("~"), ".ssh")
CONFIG = os.path.join(SSH_DIR, "config")

_ALIAS = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def existing_aliases(path=None):
    from remote.discover import ssh_config_hosts
    return ssh_config_hosts(path or CONFIG)


def block(alias, hostname, user, port=22, identity=None):
    lines = ["Host {}".format(alias), "    HostName {}".format(hostname), "    User {}".format(user)]
    if int(port) != 22:
        lines.append("    Port {}".format(int(port)))
    if identity:
        lines += ["    IdentityFile {}".format(identity), "    IdentitiesOnly yes"]
    lines.append("    ServerAliveInterval 60")
    return "\n".join(lines) + "\n"


def default_identity(alias):
    return "~/.ssh/id_ed25519_{}".format(alias)


def validate(alias, hostname, user, port):
    problems = []
    if not _ALIAS.match(alias or ""):
        problems.append("the alias must be one word (letters, digits, . _ -)")
    if not hostname or any(c.isspace() for c in hostname):
        problems.append("the address is empty or has spaces")
    if not user or any(c.isspace() for c in user) or "@" in user:
        problems.append("the user name is empty or invalid")
    try:
        if not 1 <= int(port) <= 65535:
            problems.append("the port must be 1-65535")
    except (TypeError, ValueError):
        problems.append("the port must be a number")
    return problems


def add_host(alias, hostname, user, port=22, identity=None, path=None):
    """ Appends the Host block; returns (path, backup). Raises ValueError for
        invalid fields or an alias that already exists. """
    path = path or CONFIG
    problems = validate(alias, hostname, user, port)
    if problems:
        raise ValueError("; ".join(problems))
    if alias in existing_aliases(path):
        raise ValueError("'{}' is already a Host of {}: choose another alias, or edit that file.".format(alias, path))
    os.makedirs(os.path.dirname(path), mode=0o700, exist_ok=True)
    backup = None
    old = ""
    if os.path.exists(path):
        backup = "{}.easyhybrid-{}.bak".format(path, time.strftime("%Y%m%d-%H%M%S"))
        shutil.copy2(path, backup)
        with open(path) as f:
            old = f.read()
    text = old + ("" if not old or old.endswith("\n") else "\n") + ("\n" if old.strip() else "")
    text += "# . added by EasyHybrid (Remote Machines > Add SSH Host)\n" + block(alias, hostname, user, port, identity)
    with open(path, "w") as f:
        f.write(text)
    os.chmod(path, 0o600)
    return path, backup


def key_exists(identity):
    return os.path.exists(os.path.expanduser(identity)) if identity else False


def setup_commands(alias, hostname, user, port=22, identity=None):
    """ What the user types in a terminal: create the key (if missing) and
        authorize it with the account's password (the only time it is asked). """
    identity = identity or default_identity(alias)
    cmds = []
    if not key_exists(identity):
        cmds.append("ssh-keygen -t ed25519 -f {}".format(identity))
    port_opt = " -p {}".format(int(port)) if int(port) != 22 else ""
    cmds.append("ssh-copy-id -i {}.pub{} {}@{}".format(identity, port_opt, user, hostname))
    cmds.append("ssh {} hostname".format(alias))
    return cmds
