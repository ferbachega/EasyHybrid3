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
# Parse bash environment file
# ---------------------------------------------------------------------

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
    """Check if required pDynamo modules are available."""

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

        if ask_yes_no("Would you like the installer to locate and load it for you?"):

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
# Check external libraries
# ---------------------------------------------------------------------

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


def check_external_libraries(auto_install=False):
    """ Checks required external Python libraries (see PYTHON_LIBRARIES
    above), PLUS the GTK3 system typelib (see check_gtk3_typelib()) --
    which is not a Python package at all and needs its own, separate
    fix/message.

    auto_install: if True (or if the user says yes when prompted, in
    interactive mode), attempts "pip install <package>" for whatever
    Python packages are missing. Never attempts to install the GTK3
    system package automatically -- that needs sudo/apt, a much bigger
    thing to do without asking explicitly every time, so it only ever
    prints the exact command to run. """

    print(c_header("\nChecking required Python packages...\n"))

    missing_pip_names = []

    for import_name, pip_name in PYTHON_LIBRARIES.items():
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

    if not missing_pip_names and gtk3_ok:
        print(c_ok("\n✓ All required Python packages and GTK 3 dependencies were found."))
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
        print(c_warn(
            "\nThe GTK 3.0 typelib is a SYSTEM package, not something pip can "
            "install. On Debian/Ubuntu, run:\n"
        ) + c_info("    sudo apt-get install gir1.2-gtk-3.0\n") +
            "(the exact package name may differ on other distributions)."
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

    if ask_yes_no("Create a desktop shortcut?"):

        desktop_file = Path(desktop) / "easyhybrid.desktop"

        with open(desktop_file, "w") as f:
            f.write(desktop_entry)
        os.chmod(desktop_file, 0o755)
        print(c_ok("✓ Desktop shortcut created successfully."))

    if ask_yes_no("Create an application menu entry?"):

        app_dir = Path.home() / ".local/share/applications"

        app_dir.mkdir(parents=True, exist_ok=True)

        menu_file = app_dir / "easyhybrid.desktop"

        with open(menu_file, "w") as f:
            f.write(desktop_entry)

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

This installer will:

  1. Check the required Python packages
  2. Locate and verify your pDynamo3 installation
  3. Build the VISMOL graphics engine
  4. Optionally create desktop and application-menu shortcuts

Run with --credits at any time to see the full citation details.\
''')

    if is_wsl():
        print(c_warn("\n" + WSL_NOTE_TEXT))

    # [EN] REORDERED (previously: install_vismol() ran FIRST, then
    # dependencies were checked afterwards). install_vismol() compiles
    # Cython extensions that need numpy/Cython/PyOpenGL already present
    # -- checking dependencies first means a missing package now produces
    # one clear "NOT FOUND" message up front, instead of the build step
    # failing with a confusing compiler/import error and the REAL cause
    # only becoming obvious several steps later.
    print(c_header("\n[Step 1/3] Python and system dependencies"))
    if not check_external_libraries(auto_install=NON_INTERACTIVE):
        print(c_error("\nInstallation cannot continue until the dependencies listed "
                       "above are resolved. Please address them and re-run this installer."))
        return False

    print(c_header("\n[Step 2/3] pDynamo3 installation"))
    check_pdynamo()

    print(c_header("\n[Step 3/3] VISMOL graphics engine"))
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
