#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: Installer
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      Checks and installs the dependencies EasyHybrid needs (pDynamo3,
#      Python libraries, and the VISMOL graphics engine's compiled
#      extensions), and optionally creates a desktop/application-menu
#      shortcut.
#


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EasyHybrid Installer
Python interface for QC/MM and molecular simulations using pDynamo3

Copyright 2022-2026 Fernando Bachega
License: GPL-3.0-or-later
"""

import os
import re
import sys
import shutil
import subprocess
import importlib
from pathlib import Path

from _version import EASYHYBRID_VERSION

# ---------------------------------------------------------------------
# Terminal colors
# ---------------------------------------------------------------------

# [EN] Plain ANSI escape codes rather than a dependency (colorama/rich) --
# this script has no other third-party requirements before dependencies
# are even checked, and adding one just for color would be circular.
# Disabled automatically (COLOR_ENABLED = False, every color_ helper then
# a no-op) in two cases where raw escape codes would just show up as
# garbage instead of color: output isn't an actual terminal (piped to a
# file/another program -- checked via sys.stdout.isatty()), or the user
# has NO_COLOR set (https://no-color.org, a deliberate opt-out
# convention this script honors rather than overriding).
class _Ansi:
    RESET   = "\033[0m"
    BOLD    = "\033[1m"
    RED     = "\033[31m"
    GREEN   = "\033[32m"
    YELLOW  = "\033[33m"
    BLUE    = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN    = "\033[36m"


COLOR_ENABLED = sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


def _colorize(text, *codes):
    if not COLOR_ENABLED:
        return text
    return "".join(codes) + str(text) + _Ansi.RESET


def c_header(text):
    """ Section headers / step markers, e.g. '[Step 1/3] ...'. """
    return _colorize(text, _Ansi.BOLD, _Ansi.CYAN)


def c_ok(text):
    """ Success states -- 'OK', checkmarks, completed actions. """
    return _colorize(text, _Ansi.GREEN)


def c_warn(text):
    """ Non-fatal warnings -- e.g. a missing optional/system package. """
    return _colorize(text, _Ansi.YELLOW)


def c_error(text):
    """ Failures -- 'NOT FOUND', 'ERROR: ...'. """
    return _colorize(text, _Ansi.RED)


def c_bold(text):
    """ Emphasis without implying success/failure (titles, prompts). """
    return _colorize(text, _Ansi.BOLD)


def c_info(text):
    """ Informational text that isn't a plain status (paths, commands). """
    return _colorize(text, _Ansi.BLUE)


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

EASYHYBRID_HOME = Path(__file__).resolve().parent
VISMOL_PATH = EASYHYBRID_HOME / "src" / "graphics_engine"
PDYNAMO3_EH_PATH = EASYHYBRID_HOME / "src" / "pdynamo3-eh"
DYFF_1_0EH_PATH = EASYHYBRID_HOME / "src" / "dyff-1.0eh"
DYFFDOCK_PATH = EASYHYBRID_HOME / "src" / "dyffdock"

# Every git submodule, in the order they should be checked/built. Each entry
# names the module (for the init step) and, when it needs compiling, the
# path to its own install script.
SUBMODULES = {
    "src/graphics_engine": VISMOL_PATH,
    "src/pdynamo3-eh":     PDYNAMO3_EH_PATH,
    "src/dyff-1.0eh":      DYFF_1_0EH_PATH,
    "src/dyffdock":        DYFFDOCK_PATH,
}


# ---------------------------------------------------------------------
# Credits and citation
# ---------------------------------------------------------------------

# [EN] Reference paper, confirmed via web search rather than recalled
# from memory (a citation is exactly the kind of detail worth getting
# right rather than approximately right): Bachega et al., J. Chem. Inf.
# Model. 2026, 66, 3, 1286-1292, DOI 10.1021/acs.jcim.5c02047. The
# underlying pDynamo3 library has its own separate citation (Field,
# J. Chem. Inf. Model. 2022, 62, 23, 5849-5854), included as well since
# EasyHybrid is built entirely on top of it.
CREDITS_TEXT = """\
Development team
-----------------
Jose Fernando R. Bachega, Gustavo Hagen, Carlos Sequeiros-Borja,
Kai Nikklas, Jorge Chahine, Luis Fernando M. S. Timmers, and
Martin J. Field.

How to cite
-----------
If EasyHybrid contributes to your published work, please cite:

  J. F. R. Bachega, G. Hagen, C. Sequeiros-Borja, K. Nikklas,
  J. Chahine, L. F. M. S. Timmers, M. J. Field.
  "EasyHybrid: An Interactive Graphical Environment for Quantum,
  Classical and Hybrid Simulations with pDynamo3."
  J. Chem. Inf. Model. 2026, 66 (3), 1286-1292.
  DOI: 10.1021/acs.jcim.5c02047

EasyHybrid is built entirely on the pDynamo3 simulation library;
please also cite:

  M. J. Field. "pDynamo3: Molecular Modeling and Simulation Program."
  J. Chem. Inf. Model. 2022, 62 (23), 5849-5854.
  DOI: 10.1021/acs.jcim.2c01239
"""


def is_wsl():
    """ [EN] Detects whether this installer is running inside Windows
    Subsystem for Linux, so main() can print a short, WSL-specific note
    (see WSL_NOTE_TEXT below) -- this environment has a couple of real
    gotchas beyond a native Linux install (needing WSLg/WSL2 specifically
    for the GTK3+OpenGL window to display at all, and a serious
    performance penalty for reading/writing project files through the
    /mnt/c/... Windows filesystem bridge instead of the Linux one) that
    are worth flagging proactively rather than letting the user discover
    them via a confusing error or unexplained slowdown. Checks
    /proc/version for "microsoft"/"wsl", the standard, well-established
    way to detect WSL specifically (as opposed to any other Linux). """
    try:
        with open("/proc/version") as f:
            version_info = f.read().lower()
        return "microsoft" in version_info or "wsl" in version_info
    except FileNotFoundError:
        return False


def is_macos():
    """ [EN] True on macOS. Used to gate the extra checks below that only
    apply there: GTK3's Quartz backend cannot realize a working
    Gtk.GLArea at all (see MACOS_NATIVE_RENDERING_FIX.md), so on this
    platform EasyHybrid renders through a hidden GLFW window instead
    (src/graphics_engine/.../vismol_gtkwidget.py) -- which needs both the
    "glfw" pip package AND the native GLFW library (a SEPARATE, OS-level
    dependency pip cannot install, exactly like the GTK3 typelib case
    below). Neither is needed, or checked, on Linux/Windows. """
    return sys.platform == "darwin"


WSL_NOTE_TEXT = """\
Windows Subsystem for Linux (WSL) detected.

EasyHybrid needs a working GTK3 + OpenGL display, which requires
WSL2 with WSLg (WSL's built-in graphical support, included by
default on current Windows 10/11 installations reachable via
'wsl --update' from PowerShell). WSL1 is not sufficient.

Two practical recommendations:
  - Keep the EasyHybrid checkout and any structures/trajectories you
    work with inside the Linux filesystem (e.g. under ~/), not under
    /mnt/c/... -- accessing Windows-side files through that bridge is
    considerably slower than the native Linux filesystem.
  - The GTK 3.0 typelib (gir1.2-gtk-3.0) below is a common gap on a
    freshly-installed WSL distribution; the dependency check further
    below will confirm whether it is already present.
"""


MACOS_NOTE_TEXT = """\
macOS detected.

GTK3's native (Quartz) backend cannot realize a working OpenGL
GLArea at all -- doing so corrupts window compositing for the whole
app (see MACOS_NATIVE_RENDERING_FIX.md for the full story). EasyHybrid
works around this by rendering the 3D view through a hidden GLFW
window instead of GTK's own GL widget, which needs an extra
dependency beyond every other platform:

  - The "glfw" pip package (already listed below, installed the same
    way as everything else).
  - The native GLFW library itself -- a SEPARATE, OS-level dependency
    that pip cannot install (exactly like the GTK 3.0 typelib case
    below, just the macOS equivalent of it). Install it via conda/
    mamba, e.g.:
        mamba install -c conda-forge glfw
    (Homebrew's "glfw" formula also works, if you already use brew
    instead of conda/mamba.)

The dependency check further below will confirm whether both are
already present.
"""


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

# [EN] Set by main() from the --yes/-y command-line flag. Lets the whole
# installer run unattended (e.g. inside a Docker build or a CI pipeline,
# neither of which has anyone at a keyboard to answer input() prompts) --
# every ask_yes_no() call below just returns True immediately instead of
# blocking on stdin. Previously there was no such option at all.
NON_INTERACTIVE = False


def ask_yes_no(question):
    """Prompt user for yes/no answer (auto-accepts if NON_INTERACTIVE)."""
    if NON_INTERACTIVE:
        print(f"{question} (y/n): y   [--yes]")
        return True
    answer = input(f"{question} (y/n): ").strip().lower()
    return answer in ("y", "yes")


# ---------------------------------------------------------------------
# Git submodules (graphics_engine/vismol, pdynamo3-eh, dyff-1.0eh, dyffdock)
# ---------------------------------------------------------------------

def init_submodules():
    """git submodule update --init --recursive, then report which of the 4
    expected submodule directories actually ended up non-empty.

    A plain `git clone` of EasyHybrid (without --recursive) leaves every
    submodule path present but EMPTY -- an easy trap: later steps that
    build/import from them fail with confusing "not found" errors that look
    like a missing dependency rather than a forgotten clone step. Returns
    True if every submodule is present (NOT that it builds/imports -- the
    later per-submodule steps check that), False if git itself is missing
    or the command failed.
    """
    print(c_header("\nInitializing git submodules (graphics_engine, pdynamo3-eh, dyff-1.0eh, dyffdock)...\n"))

    git = shutil.which("git")
    if git is None:
        print(c_error("ERROR: git was not found. It is needed to fetch the submodules."))
        return False

    if not (EASYHYBRID_HOME / ".git").exists():
        print(c_warn("This does not look like a git checkout of EasyHybrid (no .git here) -- "
                      "skipping submodule initialization. If the submodule folders below are "
                      "empty, you will need to fetch them yourself."))
    else:
        try:
            subprocess.run(["git", "submodule", "update", "--init", "--recursive"],
                            cwd=str(EASYHYBRID_HOME), check=True)
        except subprocess.CalledProcessError as e:
            print(c_error(f"ERROR: 'git submodule update --init --recursive' failed ({e})."))
            return False

    ok = True
    for name, path in SUBMODULES.items():
        if path.is_dir() and any(path.iterdir()):
            print(f"{name:20s} : {c_ok('present')} -> {path}")
        else:
            print(f"{name:20s} : {c_error('MISSING or empty')} -> {path}")
            ok = False
    return ok


# ---------------------------------------------------------------------
# Install VISMOL
# ---------------------------------------------------------------------

def install_vismol():
    """Install the VISMOL graphics engine, and verify the build actually
    produced usable compiled extensions (previously: a failed build only
    ever printed one message and main() carried on regardless -- the
    real failure only surfaced much later, confusingly, the first time
    the app tried to import vismol.utils.c_distances or similar).

    Returns True if the compiled extensions can be imported afterwards,
    False otherwise.
    """

    print(c_header("\nBuilding the VISMOL graphics engine...\n"))

    install_script = VISMOL_PATH / "install.sh"

    if not install_script.exists():
        print(c_error("ERROR: install.sh not found."))
        return False

    try:
        subprocess.run(
            ["bash", str(install_script)],
            cwd=str(VISMOL_PATH),
            check=True
        )
    except subprocess.CalledProcessError:
        print(c_error("ERROR: VISMOL installation failed (build_ext step raised an error -- see above)."))
        return False

    # Confirm the actual compiled Cython extensions are importable, not
    # just that the build command returned a zero exit code (setup.py
    # can, in principle, "succeed" while skipping/silently failing an
    # individual extension). vismol.utils.c_distances is the one used
    # for bond auto-detection -- if this one imports, the rest built too
    # (they're all produced by the same "python3 setup.py build_ext
    # --inplace" step).
    vismol_src = str(VISMOL_PATH / "src")
    if vismol_src not in sys.path:
        sys.path.insert(0, vismol_src)
    try:
        importlib.import_module("vismol.utils.c_distances")
        print(c_ok("\nVISMOL compiled extensions: OK (vismol.utils.c_distances imports correctly)."))
        return True
    except ImportError as e:
        print(c_error(f"\nERROR: VISMOL build finished but the compiled extension did not "
                       f"import correctly ({e}). Check the build output above for the real error."))
        return False


# ---------------------------------------------------------------------
# Build pdynamo3-eh (EasyHybrid's pDynamo3 fork, bundled as a submodule)
# ---------------------------------------------------------------------

def build_pdynamo3_eh():
    """Compile src/pdynamo3-eh (its own installation/Install.py -f: the same
    C/Cython build pDynamo3 itself uses, ~2 minutes), then verify pCore,
    pBabel and pMolecule actually import from it afterwards.

    Until 2026-11 this step instead located and patched a SEPARATELY
    installed pDynamo3 (see check_pdynamo()/install_qc_modules(), still
    below for anyone with a personal paths.py pointing elsewhere). The
    bundled fork needs none of that: it already carries every patch
    EasyHybrid used to apply by hand, so this only needs to build it.

    Returns True if pCore/pBabel/pMolecule import correctly afterwards,
    False otherwise.
    """

    print(c_header("\nBuilding pdynamo3-eh (EasyHybrid's pDynamo3 fork)...\n"))

    installer = PDYNAMO3_EH_PATH / "installation" / "Install.py"
    if not installer.exists():
        print(c_error(f"ERROR: {installer} not found (submodule not initialized -- see the step above)."))
        return False

    try:
        subprocess.run(
            [sys.executable, "Install.py", "-f"],
            cwd=str(installer.parent),
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(c_error(f"ERROR: pdynamo3-eh build failed ({e}). See the compiler output above."))
        return False

    pdynamo3_eh = str(PDYNAMO3_EH_PATH)
    if pdynamo3_eh not in sys.path:
        sys.path.insert(0, pdynamo3_eh)
    try:
        # pBabel before pMolecule: a known pDynamo3 circular-import quirk
        # (pMolecule.MMModel.MMModel imports pBabel at module load time;
        # importing pMolecule first leaves pBabel's own import of pMolecule
        # symbols reaching a partially-initialized module).
        importlib.import_module("pCore")
        importlib.import_module("pBabel")
        importlib.import_module("pMolecule")
        print(c_ok("\npdynamo3-eh: OK (pCore/pBabel/pMolecule import correctly)."))
        return True
    except ImportError as e:
        print(c_error(f"\nERROR: pdynamo3-eh build finished but did not import correctly ({e}). "
                      f"Check the compiler output above for the real error."))
        return False


# ---------------------------------------------------------------------
# dyff-1.0eh and dyffdock: pure Python, no build step -- just confirm the
# submodule is actually there and importable (needs pdynamo3-eh already on
# sys.path, i.e. build_pdynamo3_eh() must run first).
# ---------------------------------------------------------------------

def check_dyff_1_0eh():
    """Sanity-import pdynamo.dyff_edits from src/dyff-1.0eh. No compiling --
    it is pure Python -- so this only catches a missing/empty submodule or
    a genuinely broken checkout, not a build failure."""
    path = str(DYFF_1_0EH_PATH)
    if path not in sys.path:
        sys.path.insert(0, path)
    try:
        # pBabel before pMolecule: see build_pdynamo3_eh()'s own comment --
        # importing this module alone, with nothing else having imported
        # pBabel yet in this process, hits the same circular-import quirk.
        importlib.import_module("pBabel")
        importlib.import_module("pdynamo.dyff_edits")
        print(f"dyff-1.0eh          : {c_ok('OK')} -> {DYFF_1_0EH_PATH}")
        return True
    except ImportError as e:
        print(f"dyff-1.0eh          : {c_error('NOT FOUND or broken')} ({e})")
        return False


def check_dyffdock():
    """Sanity-import dyffdock from src/dyffdock (one level deeper than the
    submodule root: src/dyffdock/dyffdock/)."""
    path = str(DYFFDOCK_PATH)
    if path not in sys.path:
        sys.path.insert(0, path)
    try:
        importlib.import_module("dyffdock")
        print(f"dyffdock            : {c_ok('OK')} -> {DYFFDOCK_PATH}")
        return True
    except ImportError as e:
        print(f"dyffdock            : {c_error('NOT FOUND or broken')} ({e})")
        return False


# ---------------------------------------------------------------------
# Install external QC modules (xTB, ...) into pDynamo
# ---------------------------------------------------------------------

# EasyHybrid ships QC-model modules that are NOT part of a stock pDynamo3
# install (xTB above all). They live in EasyHybrid's src/util/extras/ and, to be
# usable, must be (1) copied into pDynamo's pMolecule/QCModel/ package and
# (2) imported from that package's __init__.py. This is the manual step users
# otherwise have to remember; automating it here avoids the classic "I replaced
# the file but nothing changed" confusion (the running copy is the one inside
# pDynamo, not the one in extras).
#
# Each entry: engine label -> (source filename in extras, symbols to import).
# The import line mirrors how pDynamo imports ORCA/DFTB:
#     from .QCModelXTB import _XTBCommand, \
#                             QCModelXTB
QC_MODULES = {
    "xTB": {
        "filename": "QCModelXTB.py",
        "module":   "QCModelXTB",
        "symbols":  ["_XTBCommand", "QCModelXTB"],
    },
    # add more here later, e.g.:
    # "SPARROW": {"filename": "QCModelSPARROW.py", "module": "QCModelSPARROW",
    #             "symbols": ["_SPARROWCommand", "QCModelSPARROW"]},
}


def _locate_pdynamo_qcmodel_dir():
    """Return the Path to pDynamo's pMolecule/QCModel package, or None."""
    try:
        m = importlib.import_module("pMolecule.QCModel")
        return Path(m.__file__).parent
    except Exception:
        return None


def _qcmodel_import_block(module, symbols):
    """Build the import line(s) matching pDynamo's own __init__ style."""
    if len(symbols) == 1:
        return "from .{:s} import {:s}\n".format(module, symbols[0])
    first = symbols[0]
    rest = symbols[1:]
    # 'from .Mod import A, \'  then continuation lines '        B'
    line = "from .{:s} import {:s}".format(module, first)
    for s in rest:
        line += ", \\\n                                        {:s}".format(s)
    return line + "\n"


def _locate_pdynamo_env_script():
    """Return the Path to pDynamo's environment_bash.com, or None.

    Derived from the pMolecule package location: the env script lives at
    <pDynamo3>/installation/shellScripts/environment_bash.com, and pMolecule is
    directly under <pDynamo3>.
    """
    try:
        m = importlib.import_module("pMolecule")
        pdynamo_root = Path(m.__file__).parent.parent
        script = pdynamo_root / "installation" / "shellScripts" / "environment_bash.com"
        return script if script.exists() else None
    except Exception:
        return None


def _write_env_var_to_script(script_path, var, value):
    """Add or update 'export VAR="value"' in a bash environment script.

    If VAR is already defined in the file, its line is replaced (so re-running
    the installer updates the path instead of appending a duplicate). Otherwise
    the export is appended. Returns True on success.
    """
    script_path = Path(script_path)
    try:
        lines = script_path.read_text().splitlines()
    except Exception:
        lines = []

    new_line = 'export {:s}="{:s}"'.format(var, value)
    # match a line that sets this var (with or without 'export', ignoring
    # leading spaces), but NOT a commented-out one
    pattern = re.compile(r'^\s*(export\s+)?' + re.escape(var) + r'\s*=')

    replaced = False
    out = []
    for line in lines:
        if pattern.match(line) and not line.lstrip().startswith("#"):
            out.append(new_line)
            replaced = True
        else:
            out.append(line)

    if not replaced:
        out.append("")
        out.append("# Added by EasyHybrid installer")
        out.append(new_line)

    try:
        # keep a one-time backup the first time we touch the file
        backup = script_path.with_suffix(script_path.suffix + ".easyhybrid.bak")
        if not backup.exists():
            shutil.copy2(str(script_path), str(backup))
        script_path.write_text("\n".join(out) + "\n")
        return True
    except Exception as e:
        print(c_error("Could not write to {:s}: {:s}".format(str(script_path), str(e))))
        return False


def configure_xtb_command():
    """Offer to persist PDYNAMO3_XTBCOMMAND into pDynamo's environment script.

    Asks for authorization, asks for the xtb executable path (validating it),
    and writes 'export PDYNAMO3_XTBCOMMAND="..."' into environment_bash.com so
    the variable is set in future shells. Returns True unless the user
    authorized it but writing failed.
    """
    var = QC_ENGINE_ENV_VARS["xTB"]  # "PDYNAMO3_XTBCOMMAND"

    # If it is already set in this shell, tell the user and offer to skip.
    current = os.environ.get(var)
    if current:
        print(c_info("\n{:s} is already set -> {:s}".format(var, current)))
        if not ask_yes_no("Do you want to (re)write it into environment_bash.com anyway?"):
            return True

    if not ask_yes_no(
        "\nMay the installer save the xTB executable path into pDynamo's "
        "environment_bash.com (sets {:s} for future shells)?".format(var)):
        print(c_info("Skipped writing {:s}. You can set it manually later.".format(var)))
        return True

    script = _locate_pdynamo_env_script()
    if script is None:
        print(c_warn(
            "Could not locate pDynamo's environment_bash.com automatically. "
            "Set {:s} manually in your shell configuration.".format(var)))
        return True

    # Ask for the xtb path, validating it is a real executable. In
    # non-interactive mode there is no one to ask -- reuse the already
    # detected path (current) if there is one, otherwise skip instead of
    # calling input() with no TTY to read from (EOFError).
    if NON_INTERACTIVE:
        if current and os.path.isfile(current) and os.access(current, os.X_OK):
            xtb_path = os.path.abspath(os.path.expanduser(current))
        else:
            print(c_info("Non-interactive mode and no valid {:s} detected; skipping.".format(var)))
            return True
    else:
        xtb_path = ""
        while True:
            xtb_path = input("\nFull path to the xtb executable "
                             "(e.g. /home/user/xtb-6.7.1/bin/xtb): ").strip()
            if not xtb_path:
                print(c_info("No path given; skipping."))
                return True
            xtb_path = os.path.abspath(os.path.expanduser(xtb_path))
            if os.path.isfile(xtb_path) and os.access(xtb_path, os.X_OK):
                break
            print(c_warn("That path is not an executable file. "
                         "Please check it and try again (or leave empty to skip)."))

    if _write_env_var_to_script(script, var, xtb_path):
        print(c_ok("Saved {:s} -> {:s}".format(var, xtb_path)))
        print(c_info("in ") + str(script))
        # also set it for the current process so later checks see it
        os.environ[var] = xtb_path
        print(c_warn(
            "Open a new shell or run:  source " + str(script) + "\n"
            "for the change to take effect in your current terminal."))
        return True
    return False


def ensure_scratch_directory():
    """Create src/pdynamo3-eh/scratch if it does not exist yet.

    pdynamo3-eh's own installer (Install.py -f) already writes
    PDYNAMO3_SCRATCH=$PDYNAMO3_HOME/scratch into environment_bash.com (read at
    GUI startup by easyhybrid.py's parse_bash_env_file(), independent of the
    current shell) -- but it never creates that directory. On a fresh clone
    nothing has, so the first QC (ORCA/xTB) calculation fails trying to write
    there. Harmless to create unconditionally: it is a plain scratch folder,
    already in src/pdynamo3-eh/.gitignore.
    """
    scratch_dir = PDYNAMO3_EH_PATH / "scratch"
    if scratch_dir.is_dir():
        return True
    try:
        scratch_dir.mkdir(parents=True, exist_ok=True)
        print(c_ok(f"Created {scratch_dir} (PDYNAMO3_SCRATCH default)."))
        return True
    except Exception as e:
        print(c_warn(f"Could not create {scratch_dir}: {e}"))
        return False


def install_qc_modules():
    """Offer to install EasyHybrid's external QC modules into pDynamo.

    For each module (currently xTB) this copies the source file from
    src/util/extras/ into pDynamo's pMolecule/QCModel/ and adds an import to
    that package's __init__.py, so 'from pMolecule.QCModel import *' exposes it.

    Returns True if nothing went wrong (including the user declining); False
    only if the user asked to install but it failed.
    """

    print(c_header("\nExternal QC modules (xTB, ...) for pDynamo\n"))

    extras_dir = EASYHYBRID_HOME / "src" / "util" / "extras"
    qcmodel_dir = _locate_pdynamo_qcmodel_dir()

    if qcmodel_dir is None:
        print(c_warn(
            "Could not locate pDynamo's pMolecule/QCModel package (is pDynamo3 "
            "loaded in this shell?). Skipping automatic QC-module installation."))

        # Give the user everything they need to do it by hand. Build the exact
        # import line(s) the installer would have added, per module.
        print(c_info(
            "\nTo install the QC module(s) manually, for each module do the "
            "following two steps:\n"))
        print("  1. Copy the module file from EasyHybrid's extras into pDynamo's\n"
              "     pMolecule/QCModel/ folder. If you don't know where that is, run:\n"
              "         python3 -c \"import pMolecule.QCModel as m; import os; "
              "print(os.path.dirname(m.__file__))\"\n"
              "     then copy the file there, e.g.:")
        for engine, info in QC_MODULES.items():
            src = extras_dir / info["filename"]
            print("         # {:s}".format(engine))
            print("         cp \"{:s}\" <pMolecule/QCModel>/".format(str(src)))
        print("\n  2. Register it by adding its import to that folder's "
              "__init__.py.\n"
              "     Append these line(s) to <pMolecule/QCModel>/__init__.py:")
        for engine, info in QC_MODULES.items():
            block = _qcmodel_import_block(info["module"], info["symbols"]).rstrip("\n")
            print("         # {:s}".format(engine))
            for bl in block.splitlines():
                print("         " + bl)
        print(c_info(
            "\nAfter that, 'from pMolecule.QCModel import *' will expose the "
            "module, and EasyHybrid will be able to use the engine.\n"
            "This is exactly what the installer does automatically when it can "
            "find the pDynamo package -- so loading the pDynamo environment "
            "first (source .../environment_bash.com) and re-running this "
            "installer will also work."))
        return True

    print(c_info("pDynamo QCModel package: ") + str(qcmodel_dir))
    print(c_info("EasyHybrid extras folder: ") + str(extras_dir))

    init_path = qcmodel_dir / "__init__.py"

    ok = True
    for engine, info in QC_MODULES.items():
        src = extras_dir / info["filename"]
        dst = qcmodel_dir / info["filename"]

        # Is it already installed (file present AND imported in __init__)?
        # The import we add/look for is 'from .<module> import ...', so detect
        # exactly that (an earlier check for 'import <module>' never matched,
        # because the module name sits after 'from .', not after 'import').
        init_text = init_path.read_text() if init_path.exists() else ""
        import_marker = "from .{:s} import".format(info["module"])
        already_imported = import_marker in init_text

        if dst.exists() and already_imported:
            print(f"{engine} : {c_ok('already installed')} -> {dst}")
            continue

        if not src.exists():
            print(f"{engine} : {c_error('source not found')} ({src})")
            ok = False
            continue

        if not ask_yes_no(
            "\nInstall the {:s} QC module into pDynamo? "
            "(copies {:s} and updates QCModel/__init__.py)".format(engine, info["filename"])):
            print(c_info("Skipped {:s}.".format(engine)))
            continue

        # 1) copy the module file
        try:
            shutil.copy2(str(src), str(dst))
            print(f"{engine} : {c_ok('copied')} -> {dst}")
        except Exception as e:
            print(f"{engine} : {c_error('copy failed')} ({e})")
            ok = False
            continue

        # 2) add the import to __init__.py (if not already there)
        if not already_imported:
            try:
                block = _qcmodel_import_block(info["module"], info["symbols"])
                with open(str(init_path), "a") as f:
                    f.write("\n# Added by EasyHybrid installer: external QC module\n")
                    f.write(block)
                print(f"{engine} : {c_ok('registered in __init__.py')}")
            except Exception as e:
                print(f"{engine} : {c_error('could not update __init__.py')} ({e})")
                print(c_warn(
                    "   Add this line manually to " + str(init_path) + ":\n"
                    "       " + _qcmodel_import_block(info["module"], info["symbols"]).strip()))
                ok = False
                continue

        # 3) verify it now imports from the package
        try:
            importlib.invalidate_caches()
            mod = importlib.import_module("pMolecule.QCModel." + info["module"])
            if hasattr(mod, info["symbols"][-1]):
                print(f"{engine} : {c_ok('verified (imports correctly)')}")
            else:
                print(f"{engine} : {c_warn('installed but symbol not found on import')}")
        except Exception as e:
            print(f"{engine} : {c_warn('installed but could not verify import')} ({e})")

    return ok




def parse_bash_env_file(filepath):
    """
    Parse a bash environment file and return environment variables.
    """

    env_vars = {}

    filepath = Path(filepath)

    if not filepath.exists():
        print("Environment file not found:", filepath)
        return env_vars

    with open(filepath, "r") as f:

        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            # Remove inline comments
            line = line.split("#")[0].strip()

            # Remove export
            line = line.replace("export ", "")

            # Remove "; export VAR"
            line = re.sub(r";\s*export\s+\w+", "", line)

            if "=" not in line:
                continue

            var, value = line.split("=", 1)

            var = var.strip()
            value = value.strip()

            value = os.path.expanduser(value)
            value = os.path.expandvars(value)

            env_vars[var] = value
            os.environ[var] = value

    # Add pDynamo to Python path
    if "PDYNAMO3_HOME" in env_vars:
        sys.path.append(env_vars["PDYNAMO3_HOME"])

    return env_vars


# ---------------------------------------------------------------------
# Check pDynamo installation
# ---------------------------------------------------------------------

def check_pdynamo():
    """Check if required pDynamo modules are available.

    [EN] 2026-11: LEGACY path, only called when a personal paths.py already
    exists (a separately installed, separately maintained pDynamo3 from
    before the pdynamo3-eh submodule existed). The default path for a fresh
    install is build_pdynamo3_eh() instead, which needs none of this."""

    print(c_header("\nLocating pDynamo3 Python modules...\n"))

    modules = ["pCore", "pMolecule", "pSimulation", "pScientific"]

    missing = []
    pdynamo_path = None

    for module in modules:

        try:

            m = importlib.import_module(module)

            module_path = Path(m.__file__).parent

            print(f"{module} : {c_ok('OK')} -> {module_path}")

            if module == "pCore":
                pdynamo_path = module_path.parent

        except ImportError:

            print(f"{module} : {c_error('NOT FOUND')}")
            missing.append(module)

    if missing:

        print(c_warn("\npDynamo3 could not be located.\n"))

        msg = """\
This does not necessarily mean pDynamo3 is not installed -- if it
is, its environment variables may simply not be loaded in this shell.

To load them manually, run:

    source <pDynamo3 installation>/installation/shellScripts/environment_bash.com
"""

        print(msg)

        if NON_INTERACTIVE:
            print(c_info("Non-interactive mode; cannot ask for the pDynamo3 path. Skipping."))

        elif ask_yes_no("Would you like the installer to locate and load it for you?"):

            pdynamo_path = input(
                "\nPath to the pDynamo3 installation (e.g. /home/user/pDynamo3): "
            ).strip()

            shell_script = Path(pdynamo_path) / "installation/shellScripts/environment_bash.com"

            if shell_script.exists():
                print(c_ok("Found environment script:"), shell_script)
                parse_bash_env_file(shell_script)

            else:
                pdynamo_path = None
                print(c_error("Could not find an environment script at that location."))

    if pdynamo_path:

        print(c_ok("\nSaving pDynamo3 path for future sessions..."))

        path_file = EASYHYBRID_HOME / "paths.py"

        with open(path_file, "w") as f:
            f.write(f'PDYNAMO_HOME = "{pdynamo_path}"\n')

    return missing


# ---------------------------------------------------------------------
# Check external QC engines (ORCA / xTB)
# ---------------------------------------------------------------------

# Each external QC engine that EasyHybrid can drive is configured through a
# pDynamo environment variable that must point at the engine executable. This
# table maps a human-readable name to that variable; add a row here when a new
# engine is supported and the check below covers it automatically.
QC_ENGINE_ENV_VARS = {
    "ORCA": "PDYNAMO3_ORCACOMMAND",
    "xTB":  "PDYNAMO3_XTBCOMMAND",
}

# The shared scratch directory both engines write their temporary files to.
QC_SCRATCH_ENV_VAR = "PDYNAMO3_SCRATCH"


def check_qc_engines():
    """Check that the external QC engines (ORCA, xTB) are configured.

    This first stage verifies the *environment variables* that pDynamo uses to
    locate each engine executable, plus the shared scratch directory. It does
    NOT yet run the executables -- that deeper check can be layered on later.

    Returns a dict {engine_name: bool} indicating whether each engine's
    environment variable is set (scratch is reported but not included in the
    per-engine result).
    """

    print(c_header("\nChecking external QC engines (environment variables)...\n"))

    results = {}

    for engine, var in QC_ENGINE_ENV_VARS.items():
        value = os.environ.get(var)
        if value:
            print(f"{engine:5s} ({var}) : {c_ok('SET')} -> {value}")
            results[engine] = True
        else:
            print(f"{engine:5s} ({var}) : {c_warn('NOT SET')}")
            results[engine] = False

    # Shared scratch directory (used by every engine).
    scratch = os.environ.get(QC_SCRATCH_ENV_VAR)
    if scratch:
        print(f"scratch ({QC_SCRATCH_ENV_VAR}) : {c_ok('SET')} -> {scratch}")
    else:
        print(f"scratch ({QC_SCRATCH_ENV_VAR}) : {c_warn('NOT SET')}")

    # Friendly summary: not finding an engine is only a warning, since a given
    # user may legitimately use just one of them (or neither).
    configured = [e for e, ok in results.items() if ok]
    if configured:
        print(c_info("\nConfigured QC engine(s): ") + ", ".join(configured))
    else:
        print(c_warn(
            "\nNo external QC engine environment variable is set. If you plan "
            "to run ORCA or xTB calculations, set the variable(s) above to the "
            "engine executable, e.g.:\n"
            "    export PDYNAMO3_XTBCOMMAND=/path/to/xtb\n"
            "    export PDYNAMO3_ORCACOMMAND=/path/to/orca"))

    if not scratch:
        print(c_warn(
            "\n" + QC_SCRATCH_ENV_VAR + " is not set. QC calculations write "
            "temporary files there; set it to a writable directory, e.g.:\n"
            "    export " + QC_SCRATCH_ENV_VAR + "=$PDYNAMO3_HOME/scratch"))

    return results




# [EN] Maps the Python import name (what importlib.import_module() needs)
# to the actual pip package name (what "pip install ..." needs) -- these
# differ for several of these (e.g. import "OpenGL" but "pip install
# PyOpenGL"). Built from direct experience setting this project up from
# scratch: this list was previously missing "gi" (PyGObject, the GTK
# bindings -- the app cannot even open a window without it) and "Cython"
# (required to build the VISMOL extensions below, but was only ever
# checked AFTER attempting that build -- see the reordered main() at the
# bottom of this file). "logging" was removed from the old list: it is
# part of the Python standard library, always present, and listing it
# alongside genuinely-external packages was misleading.
PYTHON_LIBRARIES = {
    "numpy":    "numpy",
    "OpenGL":   "PyOpenGL",
    "freetype": "freetype-py",
    "cairo":    "pycairo",
    "gi":       "PyGObject",
    "Cython":   "Cython",
    # [EN] needed to BUILD src/pdynamo3-eh (its installation/Install.py
    # reads YAML parameter files) and, at runtime, by pdynamo/
    # ligand_chemistry.py (LigandParameters.write_pdb_component).
    "yaml":     "PyYAML",
    # [EN] vismol_gtkwidget.py hard-imports this at module level (image
    # filtering/compositing for the 3D view) -- was missing from this
    # dict and from requirements.txt, so a fresh install this checker
    # reported clean still failed at first launch with
    # "ModuleNotFoundError: No module named 'PIL'".
    "PIL":      "Pillow",
    # [EN] Added for the Process Manager's "Abort" feature -- needs to
    # find and signal every DESCENDANT of a running job's process (its
    # own multiprocessing.Pool workers, or external QM programs it
    # shelled out to), not just the one direct child multiprocessing.
    # Process.terminate() already knows about. Cross-platform by
    # construction (implemented per-OS internally by psutil itself --
    # /proc on Linux, sysctl/libproc on macOS, a different mechanism
    # again on Windows), unlike an earlier version of this feature that
    # read /proc directly and consequently never worked on macOS at all.
    "psutil":   "psutil",
    # [EN] pdynamo/residue_library.py's match_residues() hard-imports this
    # (graph isomorphism against the DYFF residue-library templates --
    # Prepare DYFF System / "DYFF + residue library", whole-protein prep).
    # Was missing from this dict and from requirements.txt.
    "networkx": "networkx",
}

# [EN] macOS-only, checked in ADDITION to PYTHON_LIBRARIES above (see
# is_macos()/MACOS_NOTE_TEXT) -- never required, and never checked, on
# Linux/Windows, which keep using GTK's own Gtk.GLArea directly and have
# no use for GLFW at all.
MACOS_PYTHON_LIBRARIES = {
    "glfw": "glfw",
}


def check_gtk3_typelib():
    """ [EN] "import gi" succeeding is NOT enough to confirm the GUI can
    actually open -- gi (PyGObject) is a pip package, but the GTK 3.0
    "typelib" it introspects at runtime (gir1.2-gtk-3.0) is a SEPARATE,
    OPERATING-SYSTEM package that pip cannot install at all. Found this
    distinction the hard way, setting this project up in an environment
    that had gi installed via pip but not the system GTK3 typelib:
    "import gi" succeeds, then "gi.require_version('Gtk', '3.0')" raises
    "ValueError: Namespace Gtk not available" -- a confusing error with
    no obvious link back to "install a system package", unless you
    already know to look for it. Checked here, separately from the
    plain-import checks above/below, so the installer can give the exact
    right advice (apt, not pip) instead of a generic "gi: NOT FOUND" that
    would send the user towards the wrong fix.

    Returns True if fully available, False otherwise (including if gi
    itself isn't installed at all). """
    try:
        import gi
        gi.require_version("Gtk", "3.0")
        from gi.repository import Gtk  # noqa: F401
        return True
    except (ImportError, ValueError):
        return False


def check_glfw_native_lib():
    """ [EN] macOS-only, same distinction as check_gtk3_typelib() above,
    just for GLFW instead of GTK3: "import glfw" (the pip package,
    already checked via MACOS_PYTHON_LIBRARIES) succeeding is NOT enough
    -- it only provides Python ctypes bindings, which still need to
    locate and load the actual native GLFW shared library at runtime,
    a SEPARATE, OS-level dependency pip cannot install (see
    MACOS_NOTE_TEXT). "import glfw" succeeds either way; only
    glfw.init() actually touches the native library, so that's what's
    called here, then immediately torn down again -- this function's
    only job is the check itself, not to leave GLFW initialized.

    Returns True if the native library loads and initializes
    correctly, False otherwise (including if the "glfw" pip package
    itself isn't installed at all). """
    try:
        import glfw
        if not glfw.init():
            return False
        glfw.terminate()
        return True
    except (ImportError, AttributeError, OSError):
        return False


def check_xcode_clt():
    """ [EN] macOS-only. install_vismol() (Step 3) compiles VISMOL's
    Cython extensions (see setup.py/install.sh: "python3 setup.py
    build_ext --inplace"), which needs a working C compiler -- on
    macOS that's clang, shipped as part of the Xcode Command Line
    Tools, NOT installed by default on a fresh machine. Checked here,
    alongside the other prerequisites, for the same reason the
    dependency-check step was reordered to run before install_vismol()
    in the first place (see the comment in main()): surface one clear
    "NOT FOUND" message up front instead of letting the Cython build
    fail with a confusing "command not found: clang"/linker error deep
    in Step 3.

    "xcode-select -p" is the standard, documented way to check this: it
    prints the active developer directory and exits 0 if the Command
    Line Tools (or full Xcode) are installed, or exits with an error
    and no output if they aren't -- no compilation attempted, just asks
    Xcode's own tooling directly.

    Returns True if found, False otherwise. Never attempts to install
    them itself: "xcode-select --install" pops up Apple's own GUI
    installer and requires accepting a license there, not something
    scriptable/unattended the way a pip or conda install is. """
    try:
        result = subprocess.run(["xcode-select", "-p"], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def check_external_libraries(auto_install=False):
    """ Checks required external Python libraries (see PYTHON_LIBRARIES
    above, plus MACOS_PYTHON_LIBRARIES on macOS), PLUS the GTK3 system
    typelib (see check_gtk3_typelib()) and, on macOS, the native GLFW
    library (see check_glfw_native_lib()) and the Xcode Command Line
    Tools (see check_xcode_clt()) -- none of which are Python packages
    at all and each need their own, separate fix/message.

    auto_install: if True (or if the user says yes when prompted, in
    interactive mode), attempts "pip install <package>" for whatever
    Python packages are missing. Never attempts to install GTK3, GLFW
    or the Xcode Command Line Tools automatically -- those need sudo/
    apt, conda/brew, or Apple's own installer dialog respectively, a
    much bigger thing to do without asking explicitly every time, so it
    only ever prints the exact command to run. """

    print(c_header("\nChecking required Python packages...\n"))

    missing_pip_names = []

    # [EN] glfw is only ever needed on macOS (see MACOS_PYTHON_LIBRARIES/
    # is_macos()) -- merged into the same dict/loop as everything else so
    # it gets the exact same "NOT FOUND" reporting and pip-install flow
    # below, without duplicating either.
    libraries_to_check = dict(PYTHON_LIBRARIES)
    if is_macos():
        libraries_to_check.update(MACOS_PYTHON_LIBRARIES)

    for import_name, pip_name in libraries_to_check.items():
        try:
            importlib.import_module(import_name)
            print(f"{import_name} : {c_ok('OK')}")
        except ImportError:
            print(f"{import_name} : {c_error('NOT FOUND')}  (pip package: {pip_name})")
            missing_pip_names.append(pip_name)

    gtk3_ok = check_gtk3_typelib()
    if gtk3_ok:
        print(f"GTK 3.0 typelib (gir1.2-gtk-3.0) : {c_ok('OK')}")
    else:
        print(f"GTK 3.0 typelib (gir1.2-gtk-3.0) : {c_error('NOT FOUND')}")

    # [EN] Same "pip package imports fine, but the actual OS-level
    # library it wraps is missing" distinction as check_gtk3_typelib()
    # above, just for GLFW -- see check_glfw_native_lib(). Not applicable
    # (treated as satisfied) outside macOS.
    glfw_native_ok = True
    if is_macos():
        glfw_native_ok = check_glfw_native_lib()
        if glfw_native_ok:
            print(f"GLFW native library : {c_ok('OK')}")
        else:
            print(f"GLFW native library : {c_error('NOT FOUND')}")

    # [EN] install_vismol() (Step 3) needs a working C compiler (clang,
    # via Xcode CLT on macOS) to build VISMOL's Cython extensions -- see
    # check_xcode_clt(). Not applicable (treated as satisfied) outside
    # macOS, where a compiler is either already present (most Linux
    # dev setups) or its absence surfaces as its own clear error from
    # the build step, unrelated to anything macOS-specific here.
    xcode_clt_ok = True
    if is_macos():
        xcode_clt_ok = check_xcode_clt()
        if xcode_clt_ok:
            print(f"Xcode Command Line Tools : {c_ok('OK')}")
        else:
            print(f"Xcode Command Line Tools : {c_error('NOT FOUND')}")

    if not missing_pip_names and gtk3_ok and glfw_native_ok and xcode_clt_ok:
        print(c_ok("\n✓ All required Python packages and system dependencies were found."))
        return True

    if missing_pip_names:
        print(c_warn("\nThe following Python packages are required but were not found:"))
        for p in missing_pip_names:
            print(" -", p)
        install_cmd = ["pip", "install"] + missing_pip_names + ["--user"]
        print("\nYou can install them with:\n    " + c_info(" ".join(install_cmd)))

        do_install = auto_install or ask_yes_no(
            "\nWould you like the installer to run that pip command for you?"
        )
        if do_install:
            # [EN] Captures output (instead of letting it stream straight
            # to the terminal, as before) so the "externally-managed-
            # environment" case below can actually be detected -- found
            # by testing this exact code path, not by reading pip's
            # documentation: several modern Debian/Ubuntu systems (PEP
            # 668) refuse a plain "pip install" outside a venv, failing
            # with that specific message and suggesting either a venv or
            # --break-system-packages. Without handling it, the installer
            # would just report a bare, unexplained "pip install failed"
            # here, identical to any other pip failure.
            result = subprocess.run(install_cmd, capture_output=True, text=True)
            if result.returncode == 0:
                print(result.stdout)
                importlib.invalidate_caches()  # newly-installed package otherwise invisible to this same running process
                print(c_ok("\npip install finished -- re-checking...\n"))
                return check_external_libraries(auto_install=auto_install)

            if "externally-managed-environment" in (result.stderr or ""):
                print(c_warn(
                    "\nThis Python installation is 'externally managed' (PEP 668) "
                    "and refuses a plain pip install outside a virtual environment."
                ))
                retry_cmd = install_cmd + ["--break-system-packages"]
                if auto_install or ask_yes_no(
                    "Retry with --break-system-packages? (or answer 'n' and set up "
                    "a venv yourself instead)"
                ):
                    retry = subprocess.run(retry_cmd, capture_output=True, text=True)
                    if retry.returncode == 0:
                        print(retry.stdout)
                        importlib.invalidate_caches()  # newly-installed package otherwise invisible to this same running process
                        print(c_ok("\npip install finished -- re-checking...\n"))
                        return check_external_libraries(auto_install=auto_install)
                    print(retry.stderr)
                    print(c_error("\nERROR: pip install (with --break-system-packages) still failed."))
                    return False
                print(
                    "\nTo install manually in a virtual environment instead:\n"
                    "    python3 -m venv ~/.easyhybrid_venv\n"
                    "    source ~/.easyhybrid_venv/bin/activate\n"
                    "    " + c_info(" ".join(install_cmd))
                )
                return False

            print(result.stderr)
            print(c_error("\nERROR: pip install failed. Please install the packages manually."))
            return False

    if not gtk3_ok:
        if is_macos():
            print(c_warn(
                "\nGTK 3 (with its Python/gi bindings) was not found. This is a "
                "SYSTEM/conda package, not something pip can install. Install it "
                "via:\n"
            ) + c_info("    mamba install -c conda-forge pygobject gtk3\n") +
                '(Homebrew\'s "pygobject3" + "gtk+3" formulas also work, if you '
                "use brew instead of conda/mamba)."
            )
        else:
            print(c_warn(
                "\nThe GTK 3.0 typelib is a SYSTEM package, not something pip can "
                "install. On Debian/Ubuntu, run:\n"
            ) + c_info("    sudo apt-get install gir1.2-gtk-3.0\n") +
                "(the exact package name may differ on other distributions)."
            )

    if not glfw_native_ok:
        print(c_warn(
            "\nThe native GLFW library is a SYSTEM/conda package, not something "
            "pip can install (the 'glfw' Python package above only provides "
            "bindings to it -- see MACOS_NOTE_TEXT). Install it via:\n"
        ) + c_info("    mamba install -c conda-forge glfw\n") +
            '(Homebrew\'s "glfw" formula also works, if you use brew instead '
            "of conda/mamba)."
        )

    if not xcode_clt_ok:
        print(c_warn(
            "\nThe Xcode Command Line Tools (needed to compile VISMOL's Cython "
            "extensions in Step 3) were not found. Install them with:\n"
        ) + c_info("    xcode-select --install\n") +
            "(pops up Apple's own installer dialog -- accept the license there, "
            "then re-run this installer)."
        )

    return False



# ---------------------------------------------------------------------
# Create desktop icon
# ---------------------------------------------------------------------

def create_desktop_icon():

    try:

        desktop = subprocess.check_output(
            ["xdg-user-dir", "DESKTOP"]
        ).decode().strip()

    except Exception:

        desktop = str(Path.home() / "Desktop")

    print("\nDetected desktop directory:", desktop)

    exec_path = EASYHYBRID_HOME / "easyhybrid.py"

    icon_path = (
        EASYHYBRID_HOME
        / "src"
        / "gui"
        / "icons"
        / "easyhybrid_solo2_100x100.png"
    )

    desktop_entry = f"""[Desktop Entry]
Version=1.0
Type=Application
Name=EasyHybrid
Comment=EasyHybrid Molecular Modeling Environment
Exec={exec_path}
Icon={icon_path}
Terminal=true
Categories=Science;Education;
"""

    def _write_shortcut(path):
        """Write desktop_entry to path, backing up a pre-existing shortcut that
        pointed somewhere else first (e.g. a different EasyHybrid checkout) --
        found the hard way: running this installer against a throwaway test
        clone silently overwrote the user's real desktop/menu shortcut with
        one pointing into that (temporary) clone."""
        if path.exists():
            try:
                old = path.read_text()
            except Exception:
                old = ""
            if "Exec=" + str(exec_path) not in old:
                backup = path.with_suffix(path.suffix + ".bak")
                try:
                    shutil.copy2(str(path), str(backup))
                    print(c_warn(f"An existing shortcut at {path} pointed elsewhere; "
                                  f"backed it up to {backup}."))
                except Exception as e:
                    print(c_warn(f"Could not back up the existing {path}: {e}"))
        with open(path, "w") as f:
            f.write(desktop_entry)

    if ask_yes_no("Create a desktop shortcut?"):

        desktop_file = Path(desktop) / "easyhybrid.desktop"
        _write_shortcut(desktop_file)
        os.chmod(desktop_file, 0o755)
        print(c_ok("✓ Desktop shortcut created successfully."))

    if ask_yes_no("Create an application menu entry?"):

        app_dir = Path.home() / ".local/share/applications"

        app_dir.mkdir(parents=True, exist_ok=True)

        menu_file = app_dir / "easyhybrid.desktop"
        _write_shortcut(menu_file)

        print(c_ok("✓ Application menu entry created successfully."))


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    global NON_INTERACTIVE

    import argparse
    parser = argparse.ArgumentParser(
        description="EasyHybrid installer and dependency checker."
    )
    parser.add_argument(
        "-y", "--yes", action="store_true",
        help="Non-interactive: auto-accept every prompt (e.g. for Docker/CI builds).",
    )
    parser.add_argument(
        "--credits", action="store_true",
        help="Print author credits and the reference citation, then exit.",
    )
    parser.add_argument(
        "--profile", choices=["ff", "dock", "full"], default="full",
        help="What to install: 'ff' = pdynamo3-eh + dyff-1.0eh only (scripting/HPC use "
             "of the force field, no GUI); 'dock' = + dyffdock; 'full' (default) = "
             "+ the VISMOL graphics engine and the EasyHybrid GUI itself.",
    )
    args = parser.parse_args()

    if args.credits:
        print(c_bold(f"\nEasyHybrid {EASYHYBRID_VERSION}\n"))
        print(CREDITS_TEXT)
        return True

    NON_INTERACTIVE = args.yes

    print(c_bold(f"\nEasyHybrid {EASYHYBRID_VERSION} -- Installation and Configuration\n"))
    print('''\
EasyHybrid is a free, open-source graphical environment for molecular
simulations, built on the pDynamo3 simulation library and developed
by Bachega, Hagen, Sequeiros-Borja, Nikklas, Chahine, Timmers, and
Field (J. Chem. Inf. Model. 2026, 66, 3, 1286-1292).

It provides an intuitive interface for preparing, editing, visualizing,
and running molecular simulations, while preserving the flexibility of
Python scripting for advanced workflows.

This installer will (--profile {ff,dock,full}, default full):

  1. Check the required Python packages
  2. Fetch and build pdynamo3-eh (EasyHybrid's bundled pDynamo3 fork) and
     dyff-1.0eh (the DYFF-1.0eh force field) -- every profile
  3. [dock, full] Fetch dyffdock (docking)
  4. [full] Build the VISMOL graphics engine; configure QC engines (ORCA/xTB);
     optionally create desktop and application-menu shortcuts

Run with --credits at any time to see the full citation details.\
''')

    if is_wsl():
        print(c_warn("\n" + WSL_NOTE_TEXT))

    if is_macos():
        print(c_warn("\n" + MACOS_NOTE_TEXT))

    profile = args.profile
    print(c_info(f"\nProfile: {profile}"))

    # [EN] REORDERED (previously: install_vismol() ran FIRST, then
    # dependencies were checked afterwards). install_vismol() compiles
    # Cython extensions that need numpy/Cython/PyOpenGL already present
    # -- checking dependencies first means a missing package now produces
    # one clear "NOT FOUND" message up front, instead of the build step
    # failing with a confusing compiler/import error and the REAL cause
    # only becoming obvious several steps later.
    print(c_header("\n[Step 1] Python and system dependencies"))
    if not check_external_libraries(auto_install=NON_INTERACTIVE):
        print(c_error("\nInstallation cannot continue until the dependencies listed "
                       "above are resolved. Please address them and re-run this installer."))
        return False

    print(c_header("\n[Step 2] Git submodules"))
    init_submodules()

    print(c_header("\n[Step 3] pdynamo3-eh + dyff-1.0eh (the force field)"))
    # [EN] 2026-11: this used to locate and patch a SEPARATELY installed
    # pDynamo3 (check_pdynamo() + install_qc_modules(), kept below only for
    # a personal paths.py pointing elsewhere -- see check_pdynamo()'s own
    # docstring). The bundled submodule needs none of that.
    if (EASYHYBRID_HOME / "paths.py").exists():
        print(c_info("paths.py already present -- using the pDynamo3 it points to "
                      "instead of building the bundled pdynamo3-eh submodule."))
        check_pdynamo()
        install_qc_modules()
    elif not build_pdynamo3_eh():
        print(c_error("\npdynamo3-eh did not build or import correctly (see the errors "
                       "above); nothing that needs pDynamo3 will work until this is resolved."))
        return False
    else:
        ensure_scratch_directory()
    if not check_dyff_1_0eh():
        print(c_error("\ndyff-1.0eh did not import correctly (see the error above)."))
        return False

    if profile in ("dock", "full"):
        print(c_header("\n[Step 4] dyffdock (docking)"))
        check_dyffdock()

    if profile == "full":
        # External QC engines (ORCA / xTB) are optional but commonly used;
        # report their environment-variable configuration.
        check_qc_engines()

        # Offer to persist the xTB executable path into pDynamo's environment script.
        configure_xtb_command()

        print(c_header("\n[Step 5] VISMOL graphics engine"))
        if not install_vismol():
            print(c_error("\nVISMOL did not build or import correctly (see the errors above); "
                           "EasyHybrid's 3D visualization will not work until this is resolved."))
            return False

        create_desktop_icon()

    print(c_ok("\nEasyHybrid has been installed successfully.\n"))
    print(CREDITS_TEXT)
    return True


# ---------------------------------------------------------------------

if __name__ == "__main__":
    main()
