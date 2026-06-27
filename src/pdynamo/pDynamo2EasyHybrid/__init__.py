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
#  pDynamo2EasyHybrid/  --  layout apos refatoracao
#
#  src/pdynamo/
#  +-- pDynamo2EasyHybrid/        (pacote; __init__.py e a fachada, re-exporta
#      |                           pDynamoSession p/ manter o import antigo)
#      +-- session.py              pDynamoSession  (junta os 6 mixins via heranca)
#      +-- simulations_mixin.py    pSimulations
#      +-- analysis_mixin.py       pAnalysis
#      +-- representation.py       ModifyRepInVismol
#      +-- io_data.py              LoadAndSaveData
#      +-- import_trajectory.py    EasyHybridImportTrajectory
#      +-- restraints.py           Restraints
#      +-- helpers.py              Atom, generate_random_code, export_special_PDB
# -----------------------------------------------------------------------------

from pdynamo.pDynamo2EasyHybrid.session import pDynamoSession
from pdynamo.pDynamo2EasyHybrid.simulations_mixin import pSimulations
from pdynamo.pDynamo2EasyHybrid.analysis_mixin import pAnalysis
from pdynamo.pDynamo2EasyHybrid.representation import ModifyRepInVismol
from pdynamo.pDynamo2EasyHybrid.io_data import LoadAndSaveData
from pdynamo.pDynamo2EasyHybrid.import_trajectory import EasyHybridImportTrajectory
from pdynamo.pDynamo2EasyHybrid.restraints import Restraints
from pdynamo.pDynamo2EasyHybrid.helpers import Atom, generate_random_code, export_special_PDB
