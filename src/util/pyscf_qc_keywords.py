#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: PySCF method/functional/basis-set keyword lists for the QC Setup GUI
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Mirrors util/orca_qc_keywords.py's style (a short description inline,
#      the SetupPySCFWindow splits on the first space to get the keyword
#      itself) but with PySCF's own, lowercase naming convention -- these
#      strings are passed straight to PySCF (QCModelPySCF.orbitalBasis /
#      .functional / .method), not translated. A smaller, curated list than
#      ORCA's (PySCF's basis/functional parsers accept many more aliases
#      than listed here -- this covers what most QC/MM users actually
#      reach for; extend it if something commonly used is missing).
#

type_of_method = {
    0: "Hartree-Fock",
    1: "DFT",
}

# . The actual PySCF "method" string (QCModelPySCF.method, passed to
#   mol.apply(method)) is built from this + the main window's Restricted/
#   Unrestricted radio button: RHF/UHF for Hartree-Fock, RKS/UKS for DFT
#   (SetupPySCFWindow.on_combobox_pyscf_type_changed), the same pattern
#   SetupORCAWindow uses for restrited_label.

functionals = {
    0: "LDA LDA,VWN (Slater exchange + VWN5 correlation)",
    1: "PBE Perdew-Burke-Ernzerhof GGA",
    2: "BLYP Becke '88 exchange + Lee-Yang-Parr correlation",
    3: "BP86 Becke '88 exchange + Perdew '86 correlation",
    4: "B3LYP The popular 3-parameter hybrid (20% HF exchange)",
    5: "PBE0 One-parameter hybrid version of PBE (25% HF exchange)",
    6: "CAM-B3LYP Range-separated hybrid (better for charge transfer/excited states)",
    7: "wB97X-D Range-separated hybrid with empirical dispersion",
    8: "M06-2X Minnesota hybrid meta-GGA (main-group thermochemistry/kinetics)",
    9: "TPSS Meta-GGA (no exact exchange)",
    10: "TPSSh Hybrid meta-GGA version of TPSS (10% HF exchange)",
    11: "SCAN Meta-GGA (Strongly Constrained and Appropriately Normed)",
}

basis_set_group_dict = {
    0: "Pople style basis sets",
    1: "Ahlrichs def2 basis sets",
    2: "Dunning correlation-consistent basis sets",
    3: "Minimal/semi-empirical-like basis sets",
}

list_Pople_Basis = [
    "6-31G Pople 6-31G",
    "6-31G* 6-31G plus one polarisation function on non-hydrogen atoms",
    "6-31G** 6-31G plus one polarisation function on all atoms",
    "6-31+G* 6-31G* plus diffuse functions on non-hydrogen atoms",
    "6-31++G** 6-31G** plus diffuse functions on all atoms",
    "6-311G 6-311G",
    "6-311G** 6-311G plus one polarisation function on all atoms",
    "6-311+G** 6-311G** plus diffuse functions on non-hydrogen atoms",
    "6-311++G** 6-311G** plus diffuse functions on all atoms",
]

list_Def2_Basis = [
    "def2-SVP Ahlrichs def2 split-valence + polarisation",
    "def2-SVPD def2-SVP plus diffuse functions",
    "def2-TZVP Ahlrichs def2 triple-zeta valence + polarisation",
    "def2-TZVPD def2-TZVP plus diffuse functions",
    "def2-TZVPP def2-TZVP with extra polarisation functions",
    "def2-QZVP Ahlrichs def2 quadruple-zeta valence + polarisation",
]

list_Dunning_Basis = [
    "cc-pVDZ Dunning correlation-consistent double-zeta",
    "cc-pVTZ Dunning correlation-consistent triple-zeta",
    "cc-pVQZ Dunning correlation-consistent quadruple-zeta",
    "aug-cc-pVDZ cc-pVDZ plus diffuse functions",
    "aug-cc-pVTZ cc-pVTZ plus diffuse functions",
]

list_Minimal_Basis = [
    "STO-3G Minimal (Slater-type orbitals fitted to 3 Gaussians)",
    "3-21G Pople 3-21G (split-valence, no polarisation)",
]

type_basis_dict = {
    0: list_Pople_Basis,
    1: list_Def2_Basis,
    2: list_Dunning_Basis,
    3: list_Minimal_Basis,
}
