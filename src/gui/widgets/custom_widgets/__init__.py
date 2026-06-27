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
# -----------------------------------------------------------------------------
#  custom_widgets/  --  layout apos refatoracao
#
#  src/gui/widgets/
#  +-- custom_widgets/             (pacote; __init__.py e a fachada, re-exporta
#      |                            tudo, mantendo os imports antigos funcionando)
#      +-- pixmap_helpers.py        getColouredPixmap, get_colorful_square_pixel_buffer
#      +-- file_choosers.py         FolderChooserButton, FileChooser
#      +-- save_trajectory_box.py   SaveTrajectoryBox
#      +-- selection_box.py         VismolSelectionTypeBox
#      +-- trajectory_frame.py      VismolTrajectoryFrame
#      +-- comboboxes.py            SystemComboBox, CoordinatesComboBox
#      +-- reaction_coordinate.py   ReactionCoordinateBox, AdvancedReactionCoordinateBox,
#      |                            compute_sigma_a1_a3
#      +-- sequence_viewer_box.py   SequenceViewerBox
# -----------------------------------------------------------------------------

from gui.widgets.custom_widgets.pixmap_helpers import getColouredPixmap, get_colorful_square_pixel_buffer
from gui.widgets.custom_widgets.file_choosers import FolderChooserButton, FileChooser
from gui.widgets.custom_widgets.save_trajectory_box import SaveTrajectoryBox
from gui.widgets.custom_widgets.selection_box import VismolSelectionTypeBox
from gui.widgets.custom_widgets.trajectory_frame import VismolTrajectoryFrame
from gui.widgets.custom_widgets.comboboxes import SystemComboBox, CoordinatesComboBox
from gui.widgets.custom_widgets.reaction_coordinate import compute_sigma_a1_a3, AdvancedReactionCoordinateBox, ReactionCoordinateBox
from gui.widgets.custom_widgets.sequence_viewer_box import SequenceViewerBox
