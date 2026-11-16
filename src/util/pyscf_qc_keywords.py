#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: PySCF method/functional/basis-set/ECP keyword lists for the QC Setup GUI
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
#      .functional / .method / mole_kwargs['ecp']), not translated.
#
#      Every entry below was checked against a real PySCF install (2.15.0,
#      pyscf.dft.libxc.parse_xc() for functionals, pyscf.gto.basis.load()/
#      load_ecp() for basis sets/ECPs) -- 2026-11-16, see the commit this
#      file was added/expanded in for the exact check script. PySCF/libxc
#      accept many more names than listed here (this is a curated subset of
#      what QC/MM users actually reach for, not the full libxc/basis
#      library) -- SetupPySCFWindow's "Custom Functional"/"Custom Basis
#      Set"/"ECP" entries let anyone type a string beyond this list; it
#      does not need to be re-verified here to be usable, only typed
#      correctly (PySCF's own parsers are case-insensitive and ignore
#      dashes/underscores in basis names).
#

type_of_method = {
    0: "Hartree-Fock",
    1: "DFT",
}

# . The actual PySCF "method" string (QCModelPySCF.method, passed to
#   mol.apply(method)) is built from this + the main window's Restricted/
#   Unrestricted radio button: RHF/UHF for Hartree-Fock, RKS/UKS for DFT
#   (SetupPySCFWindow.on_combobox_pyscf_type_changed), the same pattern
#   SetupORCAWindow uses for restrited_label. QCModelPySCF's own module
#   docstring restricts itself to exactly these 4 mean-field methods for
#   QC/MM energies and gradients (post-HF methods such as MP2/CCSD/CASSCF
#   -- which PySCF itself supports widely -- would need to be applied by
#   hand to the converged system.qcModel.mf afterwards; this window's
#   "Type" choice cannot offer them without silently producing a system
#   the energy/gradient driver cannot actually run).

functionals = {
    0:  "LDA LDA,VWN (Slater exchange + VWN5 correlation)",
    1:  "SVWN Slater exchange + VWN5 correlation (same as LDA)",
    2:  "PBE Perdew-Burke-Ernzerhof GGA",
    3:  "PBEsol PBE revised for solids",
    4:  "BLYP Becke '88 exchange + Lee-Yang-Parr correlation",
    5:  "BP86 Becke '88 exchange + Perdew '86 correlation",
    6:  "OLYP Handy's 'optimal' exchange + Lee-Yang-Parr correlation",
    7:  "RPBE Revised PBE (improved chemisorption energies)",
    8:  "PW91 Perdew-Wang '91 GGA functional",
    9:  "revTPSS Revised TPSS meta-GGA",
    10: "TPSS Meta-GGA (Tao-Perdew-Staroverov-Scuseria, no exact exchange)",
    11: "SCAN Meta-GGA (Strongly Constrained and Appropriately Normed)",
    12: "R2SCAN Regularised, more numerically stable SCAN",
    13: "M06-L Minnesota local (no exact exchange) meta-GGA",
    14: "B3LYP The popular 3-parameter hybrid (20% HF exchange)",
    15: "PBE0 One-parameter hybrid version of PBE (25% HF exchange)",
    16: "X3LYP Extended hybrid (improves on B3LYP for H-bonding/vdW)",
    17: "O3LYP Hybrid version of OLYP",
    18: "HSE06 Screened (range-separated) hybrid, good for solids",
    19: "TPSSh Hybrid meta-GGA version of TPSS (10% HF exchange)",
    20: "M06 Minnesota hybrid meta-GGA (transition metals/main group)",
    21: "M06-2X Minnesota hybrid meta-GGA (main-group thermochemistry/kinetics)",
    22: "M05-2X Minnesota hybrid meta-GGA (an earlier M0X variant)",
    23: "MN15 Minnesota hybrid meta-GGA (broad main-group/TM accuracy)",
    24: "CAM-B3LYP Range-separated hybrid (better for charge transfer/excited states)",
    25: "wB97X Range-separated hybrid",
    26: "wB97X-D Range-separated hybrid with empirical dispersion",
    27: "wB97X-V Range-separated hybrid with VV10 non-local correlation",
    28: "wB97M-V Range-separated meta-hybrid with VV10 non-local correlation",
}

basis_set_group_dict = {
    0: "Pople style basis sets",
    1: "Ahlrichs def2 basis sets",
    2: "Dunning correlation-consistent basis sets",
    3: "Jensen polarization-consistent basis sets",
    4: "Minimal basis sets",
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
    "def2-QZVPD def2-QZVP plus diffuse functions",
    "def2-QZVPP def2-QZVP with extra polarisation functions",
]

list_Dunning_Basis = [
    "cc-pVDZ Dunning correlation-consistent double-zeta",
    "cc-pVTZ Dunning correlation-consistent triple-zeta",
    "cc-pVQZ Dunning correlation-consistent quadruple-zeta",
    "aug-cc-pVDZ cc-pVDZ plus diffuse functions (anions/weak interactions)",
    "aug-cc-pVTZ cc-pVTZ plus diffuse functions",
    "aug-cc-pVQZ cc-pVQZ plus diffuse functions",
]

list_Jensen_Basis = [
    "pc-0 Jensen polarization-consistent, minimal (DFT-optimised)",
    "pc-1 Jensen polarization-consistent, double-zeta",
    "pc-2 Jensen polarization-consistent, triple-zeta",
    "pc-3 Jensen polarization-consistent, quadruple-zeta",
    "aug-pc-1 pc-1 plus diffuse functions",
    "aug-pc-2 pc-2 plus diffuse functions",
]

list_Minimal_Basis = [
    "STO-3G Minimal (Slater-type orbitals fitted to 3 Gaussians)",
    "STO-6G Minimal (Slater-type orbitals fitted to 6 Gaussians, more accurate than STO-3G)",
    "3-21G Pople 3-21G (split-valence, no polarisation)",
]

type_basis_dict = {
    0: list_Pople_Basis,
    1: list_Def2_Basis,
    2: list_Dunning_Basis,
    3: list_Jensen_Basis,
    4: list_Minimal_Basis,
}

# . Effective core potentials (relativistic core replacement for heavy
#   elements) -- optional, wired through QCModelPySCF's mole_kwargs
#   ("ecp"). "None" (index 0) means all-electron (no ECP), PySCF's own
#   default. A matching orbital basis (e.g. def2-* with the def2-* ECP of
#   the same name, or the -PP correlation-consistent sets with their own
#   ECP) is the user's responsibility -- same as in any PySCF script.
ecp_dict = {
    0: "None (all-electron)",
    1: "def2-SVP def2 ECP (use with a def2-SVP/def2-TZVP/... orbital basis)",
    2: "def2-TZVP def2 ECP (use with a def2-TZVP/def2-QZVP/... orbital basis)",
    3: "LANL2DZ Los Alamos ECP, double-zeta valence",
    4: "LANL2TZ Los Alamos ECP, triple-zeta valence",
    5: "Stuttgart Stuttgart/Dresden small-core ECP",
    6: "CRENBL Christiansen/Ermler large-core ECP (incl. spin-orbit variants)",
    7: "SBKJC Stevens/Basch/Krauss/Jasien/Cundari ECP",
}
