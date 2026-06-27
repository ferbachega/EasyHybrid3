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


#from LogFile import LogFileWriter
# pDynamo
# -----------------------------------------------------------------------------
#  p_methods/  --  layout apos refatoracao
#
#  src/pdynamo/
#  +-- p_methods/                 (pacote; __init__.py e a fachada, re-exporta
#      |                           todas as classes p/ manter imports antigos)
#      +-- _common.py              backup_orca_files, get_hamiltonian, write_header,
#      |                           plot_data, func        (utilitarios transversais)
#      +-- energy.py               LogFile, EnergyCalculation, EnergyRefinement
#      +-- geometry_optimization.py GeometryOptimization
#      +-- molecular_dynamics.py   MolecularDynamics
#      +-- surface_scan.py         RelaxedSurfaceScan, AdvancedRelaxedSurfaceScan
#      |                           + helpers de scan/coordenada de reacao (paralelo)
#      +-- umbrella_sampling.py    UmbrellaSampling
#      |                           + funcoes _run_*_parallel_* (picklaveis p/ mp)
#      |                           + helpers _us_*
#      +-- wham.py                 WHAMAnalysis
#      +-- chain_of_states.py      ChainOfStatesOptimizePath
#      +-- normal_modes.py         NormalModes
# -----------------------------------------------------------------------------

from pdynamo.p_methods.energy import LogFile
from pdynamo.p_methods.energy import EnergyCalculation
from pdynamo.p_methods.energy import EnergyRefinement
from pdynamo.p_methods.geometry_optimization import GeometryOptimization
from pdynamo.p_methods.molecular_dynamics import MolecularDynamics
from pdynamo.p_methods.surface_scan import AdvancedRelaxedSurfaceScan
from pdynamo.p_methods.surface_scan import RelaxedSurfaceScan
from pdynamo.p_methods.umbrella_sampling import UmbrellaSampling
from pdynamo.p_methods.wham import WHAMAnalysis
from pdynamo.p_methods.chain_of_states import ChainOfStatesOptimizePath
from pdynamo.p_methods.normal_modes import NormalModes
