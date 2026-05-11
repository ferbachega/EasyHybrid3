#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
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
#      Provides functions for selecting atoms and residues in pDynamo systems
#      to facilitate QC/MM partitioning and molecular simulations.
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

EASYHYBRID_VERSION = "3.0.2"

# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

EASYHYBRID_HOME = Path(__file__).resolve().parent
VISMOL_PATH = EASYHYBRID_HOME / "src" / "graphics_engine"


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def ask_yes_no(question):
    """Prompt user for yes/no answer."""
    answer = input(f"{question} (y/n): ").strip().lower()
    return answer in ("y", "yes")


# ---------------------------------------------------------------------
# Install VISMOL
# ---------------------------------------------------------------------

def install_vismol():
    """Install the VISMOL graphics engine."""

    print("\nInstalling VISMOL Graphics Engine...\n")

    install_script = VISMOL_PATH / "install.sh"

    if not install_script.exists():
        print("ERROR: install.sh not found.")
        return

    try:
        subprocess.run(
            ["bash", str(install_script)],
            cwd=str(VISMOL_PATH),
            check=True
        )
    except subprocess.CalledProcessError:
        print("ERROR: VISMOL installation failed.")


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

    print("\n=== Checking pDynamo Python modules ===\n")

    modules = ["pCore", "pMolecule", "pSimulation", "pScientific"]

    missing = []
    pdynamo_path = None

    for module in modules:

        try:

            m = importlib.import_module(module)

            module_path = Path(m.__file__).parent

            print(f"{module} : OK -> {module_path}")

            if module == "pCore":
                pdynamo_path = module_path.parent

        except ImportError:

            print(f"{module} : NOT FOUND")
            missing.append(module)

    if missing:

        print("\n=== pDynamo not properly installed ===\n")

        msg = """
Apparently pDynamo is not properly installed.

If you have already installed it, the environment variables
may not yet be loaded.

Run the following command:

    source environment_bash.com

located in:

    ../pDynamo3/installation/shellScripts/
"""

        print(msg)

        if ask_yes_no("Would you like the installer to check this for you?"):

            pdynamo_path = input(
                "\nEnter the path to the pDynamo installation (e.g. /home/user/pDynamo3): "
            ).strip()

            shell_script = Path(pdynamo_path) / "installation/shellScripts/environment_bash.com"

            if shell_script.exists():
                print("Found environment script:", shell_script)
                parse_bash_env_file(shell_script)

            else:
                pdynamo_path = None
                print("Environment script not found.")

    if pdynamo_path:

        print("\nSaving pDynamo path file.")

        path_file = EASYHYBRID_HOME / "paths.py"

        with open(path_file, "w") as f:
            f.write(f'PDYNAMO_HOME = "{pdynamo_path}"\n')

    return missing


# ---------------------------------------------------------------------
# Check external libraries
# ---------------------------------------------------------------------

def check_external_libraries():
    """Check required external Python libraries."""

    modules = ["numpy", "OpenGL", "logging", "freetype", "cairo"]

    missing = []
    failed = False
    print("\nChecking required Python libraries:\n")

    for module in modules:

        try:

            importlib.import_module(module)

            print(f"{module} : OK")

        except ImportError:

            print(f"{module} : NOT FOUND")
            missing.append(module)
            failed = True
    
    if failed:
        print("\nMissing required libraries:\n")
        for m in missing:
            print(" -", m)

        return False
    else:
        print("\nAll required Python libraries are installed.")
        return True


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

    if ask_yes_no("Create a Desktop shortcut?"):

        desktop_file = Path(desktop) / "easyhybrid.desktop"

        with open(desktop_file, "w") as f:
            f.write(desktop_entry)
        os.chmod(desktop_file, 0o755)
        print("Desktop shortcut created.")

    if ask_yes_no("Create an application menu entry?"):

        app_dir = Path.home() / ".local/share/applications"

        app_dir.mkdir(parents=True, exist_ok=True)

        menu_file = app_dir / "easyhybrid.desktop"

        with open(menu_file, "w") as f:
            f.write(desktop_entry)

        print("Application menu entry created.")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    print(f"\nEasyHybrid Installer (version {EASYHYBRID_VERSION})\n")

    install_vismol()

    check_pdynamo()

    if not check_external_libraries():

        print("\nPlease install the missing dependencies.")
        return False
    
    create_desktop_icon()
    print("\nInstallation check completed.\n")


# ---------------------------------------------------------------------

if __name__ == "__main__":
    main()

