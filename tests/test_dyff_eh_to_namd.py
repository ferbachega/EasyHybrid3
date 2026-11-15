"""
DYFF (EasyHybrid) system -> NAMD/CHARMM files (src/util/dyff_eh_to_namd.py), 2026-10-15.

The full check (NAMD single points vs pDynamo, minimization) is
examples/dyff_eh_namd_export/convert_to_namd.py. Here: path parsing, the
one-parameter-set-per-type-tuple rule, exact torsion conversion, and the PRM
title read by namd_runner (Run NAMD window).
"""
import math
import os
import sys

import pytest

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (os.path.join(_REPO, "src"), os.path.join(_REPO, "src", "graphics_engine", "src")):
    if p not in sys.path:
        sys.path.insert(0, p)

from util import dyff_eh_to_namd as D        # noqa: E402

PKL = os.path.join(_REPO, "examples", "dyff_eh_system_preparation", "output", "valacyclovir_dyff_eh.pkl")


def test_parse_atom_path():
    assert D.parse_atom_path("A:GLN.18:N") == ("A", "GLN", 18, "N")
    assert D.parse_atom_path(":LIG.1:C12") == ("SYS", "LIG", 1, "C12")
    assert D.parse_atom_path("") == ("SYS", "UNK", 1, "X")


def test_type_names_are_short():
    assert D._type_name(0) == "D0" and len(D._type_name(36 ** 4)) <= 6


@pytest.fixture(scope="module")
def converted(tmp_path_factory):
    if not os.path.exists(PKL):
        pytest.skip("run examples/dyff_eh_system_preparation/prepare_dyff_eh.py first")
    pytest.importorskip("pBabel")
    prefix = str(tmp_path_factory.mktemp("namd") / "vala")
    system, data = D.convert(PKL, prefix)
    return system, data, prefix


def test_one_parameter_set_per_type_tuple(converted):
    system, data, prefix = converted
    for fam, items, params in (("B", data.bonds, data.bondParams), ("A", data.angles, data.angleParams),
                               ("D", data.dihedrals, data.dihedralParams)):
        for ats in items:
            assert D._canonical(tuple(data.typeOfAtom[b] for b in ats)) in params
    for ats in data.impropers:
        assert tuple(data.typeOfAtom[b] for b in ats) in data.improperK


def test_dihedrals_exact(converted):
    # . NAMD form from the PRM lines = pDynamo cosine series minus the constants
    system, data, prefix = converted
    from util.dyff_to_charmm import charmm_psi
    x = data.coords
    namd = 0.0
    for ats in data.dihedrals:
        phi = charmm_psi(*[x[b] for b in ats])
        for K, n, delta in data.dihedralParams[D._canonical(tuple(data.typeOfAtom[b] for b in ats))]:
            namd += K * (1.0 + math.cos(n * phi - math.radians(delta)))
    pdyn = D.pdynamo_energies(system)["DIHED"]
    assert abs(namd - (pdyn - data.constant_kcal)) < 1.0e-6


def test_prm_declares_namd_settings(converted):
    system, data, prefix = converted
    settings = D.namd_settings_from_prm([prefix + ".prm"])
    assert settings == {"vdwGeometricSigma": True, "1-4scaling": 0.5}
    from util import namd_runner as R
    cfg = R.build_namd_config("x", charmm_psf=prefix + ".psf", charmm_parameters=[prefix + ".prm"],
                              coordinates=prefix + ".pdb", minimize_steps=10, pbc=False)
    assert "vdwGeometricSigma               yes" in cfg and "1-4scaling                      0.5" in cfg
