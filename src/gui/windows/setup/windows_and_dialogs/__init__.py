#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Copyright 2022-2025 Fernando Bachega
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
#      to facilitate QM/MM partitioning and molecular simulations.
#
# ─────────────────────────────────────────────────────────────────────────────
#  windows_and_dialogs/  —  refactored layout
#
#  src/gui/windows/setup/
#  └── windows_and_dialogs/          (package; __init__.py is the facade —
#                                     re-exports everything below, so existing
#                                     `from ...windows_and_dialogs import X`
#                                     imports keep working unchanged)
#      ├── dialogs/
#      │   ├── simple_dialog.py       SimpleDialog, call_message_dialog
#      │   ├── text_info_windows.py   TextWindow, InfoWindow
#      │   ├── export_dialogs.py      ExportScriptDialog, ExportDataWindow
#      │   └── restraint_dialog.py    AddHarmonicRestraintDialog
#      ├── qc_setup/
#      │   ├── setup_sparrow.py       SetupSPARROWWindow
#      │   ├── setup_mopac.py         SetupMOPACWindow
#      │   ├── setup_xtb.py           SetupXTBWindow
#      │   ├── setup_dftbplus.py      SetupDFTBplusWindow
#      │   ├── setup_orca.py          SetupORCAWindow
#      │   └── qcmodel_window.py      EasyHybridSetupQCModelWindow
#      ├── system_windows/
#      │   ├── import_system.py       ImportANewSystemWindow
#      │   ├── import_trajectory.py   ImportTrajectoryWindow
#      │   ├── merge_system.py        MergeSystemWindow
#      │   ├── solvate_system.py      SolvateSystemWindow, MakeSolventBoxWindow
#      │   └── trajectory_player.py   TrajectoryPlayerWindow
#      └── selection_windows/
#          ├── selection.py          EasyHybridSelectionWindow, PDynamoSelectionWindow
#          ├── go_to_atom.py         EasyHybridGoToAtomWindow
#          ├── set_qc_atoms.py       EasyHybridDialogSetQCAtoms
#          └── prune.py              EasyHybridDialogPrune
# ─────────────────────────────────────────────────────────────────────────────

from gui.windows.setup.windows_and_dialogs.dialogs.simple_dialog import call_message_dialog, SimpleDialog
from gui.windows.setup.windows_and_dialogs.dialogs.text_info_windows import TextWindow, InfoWindow
from gui.windows.setup.windows_and_dialogs.dialogs.export_dialogs import ExportScriptDialog, ExportDataWindow
from gui.windows.setup.windows_and_dialogs.dialogs.restraint_dialog import AddHarmonicRestraintDialog
from gui.windows.setup.windows_and_dialogs.dialogs.restraint_dialog import AddPositionHarmonicRestraintDialog
from gui.windows.setup.windows_and_dialogs.qc_setup.setup_sparrow import SetupSPARROWWindow
from gui.windows.setup.windows_and_dialogs.qc_setup.setup_mopac import SetupMOPACWindow
from gui.windows.setup.windows_and_dialogs.qc_setup.setup_xtb import SetupXTBWindow
from gui.windows.setup.windows_and_dialogs.qc_setup.setup_dftbplus import SetupDFTBplusWindow
from gui.windows.setup.windows_and_dialogs.qc_setup.setup_orca import SetupORCAWindow
from gui.windows.setup.windows_and_dialogs.qc_setup.qcmodel_window import EasyHybridSetupQCModelWindow
from gui.windows.setup.windows_and_dialogs.system_windows.import_system import ImportANewSystemWindow
from gui.windows.setup.windows_and_dialogs.system_windows.import_trajectory import ImportTrajectoryWindow
from gui.windows.setup.windows_and_dialogs.system_windows.merge_system import MergeSystemWindow
from gui.windows.setup.windows_and_dialogs.system_windows.solvate_system import SolvateSystemWindow, MakeSolventBoxWindow
from gui.windows.setup.windows_and_dialogs.system_windows.trajectory_player import TrajectoryPlayerWindow
from gui.windows.setup.windows_and_dialogs.selection_windows.selection import EasyHybridSelectionWindow, PDynamoSelectionWindow
from gui.windows.setup.windows_and_dialogs.selection_windows.go_to_atom import EasyHybridGoToAtomWindow
from gui.windows.setup.windows_and_dialogs.selection_windows.set_qc_atoms import EasyHybridDialogSetQCAtoms
from gui.windows.setup.windows_and_dialogs.selection_windows.prune import EasyHybridDialogPrune
