"""
Testes do ajuste de diedros DYFF contra uma varredura xTB
(src/pdynamo/torsion_fit.py):

  1) o ajuste linear recupera coeficientes conhecidos (dados sinteticos);
  2) periodos que se cancelam (rotor simetrico) sao deixados de fora;
  3) uma ligacao em anel e' recusada;
  4) varredura real (xtb) da ligacao entre os aneis da flavona: o ajuste
     reduz o erro e a curva DYFF final e' recalculada pelo pDynamo.

Precisam do pDynamo3 (e do xtb, no teste 4).
"""
import os
import sys

import numpy as np
import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_PDYNAMO = os.environ.get("PDYNAMO3_HOME", os.path.expanduser("~/programs/pDynamo3"))
if _PDYNAMO not in sys.path:
    sys.path.insert(0, _PDYNAMO)
pytest.importorskip("pBabel")
os.environ.setdefault("PDYNAMO3_PARAMETERS", os.path.join(_PDYNAMO, "parameters"))

from pdynamo import torsion_fit as F


def _butane_like():
    """ 4 atoms in a chain, the last one rotated: dihedral 0..345 deg. """
    base = np.array([[1.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.5], [1.0, 1.0, 1.5]])
    side = {3}
    return [F.rotate_about_bond(base, 1, 2, side, a) for a in range(0, 360, 15)]


def test_fit_recovers_known_coefficients():
    geometries = _butane_like()
    quad = (0, 1, 2, 3)
    phi = np.radians([F.dihedral_angle(x, quad) for x in geometries])
    true = {1: 3.0, 2: -5.0, 3: 1.5}
    e_qm = sum(c * np.cos(n * phi) for n, c in true.items()) + 7.0
    pairs, dropped, w = F.fit_periods(geometries, e_qm, np.zeros(len(geometries)), [quad], max_period=3)
    got = {n: c for c, n in pairs}
    assert not dropped
    for n, c in true.items():
        assert abs(got[n] - c) < 1e-6


def test_symmetric_rotor_drops_cancelled_periods():
    """ three equivalent H's 120 deg apart: sum cos(n phi_d) = 0 for n = 1, 2. """
    j, k = np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.5])
    i = np.array([1.0, 0.0, -0.5])
    hs = [np.array([np.cos(t), np.sin(t), 2.0]) for t in np.radians([0, 120, 240])]
    base = np.array([i, j, k] + hs)
    geometries = [F.rotate_about_bond(base, 1, 2, {3, 4, 5}, a) for a in range(0, 360, 10)]
    quads = [(0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 2, 5)]
    phi = np.radians([[F.dihedral_angle(x, q) for q in quads] for x in geometries])
    e = 2.0 * np.cos(3 * phi).sum(axis=1)
    pairs, dropped, w = F.fit_periods(geometries, e, np.zeros(len(geometries)), quads, max_period=3)
    assert dropped == [1, 2]
    assert pairs == [(pytest.approx(2.0), 3)]


def test_ring_bond_is_refused():
    adjacency = F.adjacency_from_bonds(6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)])
    with pytest.raises(ValueError):
        F.moving_side(adjacency, 1, 2)


@pytest.fixture
def flavone(make_vismol_object):
    from pdynamo.opls import ligand
    atoms, bonds = ligand._parse_mol2(os.path.join(_TESTS_DIR, "data", "flavone.mol2"))
    elements = [a[1] for a in atoms]
    kekule, _ = ligand._kekulize(elements, bonds)
    vobj = make_vismol_object(elements, [x for (i, j, o) in kekule for x in (i, j)])
    for k, a in enumerate(atoms):
        vobj.frames[0, k] = a[2]
        vobj.atoms[k].name = "%s%d" % (elements[k], k + 1)
    vobj.manual_bonds = {(min(i, j), max(i, j)) for (i, j, o) in kekule}
    vobj.manual_bond_orders = {(min(i, j), max(i, j)): o for (i, j, o) in kekule}
    vobj.manual_aromatic_bonds = set()
    return vobj


def test_flavone_inter_ring_torsion_fit(flavone):
    from pdynamo.opls.ligand import find_xtb
    if find_xtb() is None:
        pytest.skip("xtb not found")
    from gui.windows.builder import dyff_parameters as D
    tables = D.build_tables(flavone)
    # . the dihedral row around the inter-ring bond C2-C4 (atoms 3, 5: 0-based)
    j, k = 3, 5
    row = next(r for r in tables["dihedral"] if any({q[1], q[2]} == {j, k} for q in r["instances"]))
    quad = next(q for q in row["instances"] if {q[1], q[2]} == {j, k})
    symbols, coordinates, bonds = D.scan_inputs(flavone)
    runner = F.XTBRunner(gfn=2, charge=0, multiplicity=1)
    scan = F.run_scan(symbols, coordinates, bonds, quad, step=30.0, relaxed=True, runner=runner)
    assert len(scan["geometries"]) == 12
    target = D.DYFFTorsionTarget(flavone, row, quad)
    result = F.fit_dihedral_key(scan, target.key_text, target.instances, lambda p: target.energies_for(p, scan["geometries"]),
                                current_text=target.current_text, format_value=target.format_value, max_period=4)
    print(F.report_text(result))
    assert result["rmse_after"] < result["rmse_before"]
    assert result["rmse_after"] < 2.5
    assert len(result["phi"]) == len(result["qm"]) == len(result["after"]) == 12


def test_flavone_rigid_and_two_direction_scans(flavone):
    """ rigid scan: the xTB-optimized molecule is only rotated (bond lengths
        kept); the two-direction relaxed scan keeps, at every angle, the lower
        energy; relaxing lowers the barrier. """
    from pdynamo.opls.ligand import find_xtb
    if find_xtb() is None:
        pytest.skip("xtb not found")
    from gui.windows.builder import dyff_parameters as D
    symbols, coordinates, bonds = D.scan_inputs(flavone)
    quad = (10, 5, 3, 0)
    runner = F.XTBRunner(gfn=2)
    rigid = F.run_scan(symbols, coordinates, bonds, quad, step=60.0, relaxed=False, runner=runner)
    assert len(rigid["e_qm"]) == 6
    first = rigid["geometries"][0]
    for xyz in rigid["geometries"]:
        assert abs(np.linalg.norm(xyz[3] - xyz[5]) - np.linalg.norm(first[3] - first[5])) < 1e-6
    steps = np.diff(sorted(F.wrap180(p - rigid["phi"][0]) for p in rigid["phi"]))
    assert np.allclose(steps, 60.0, atol=1e-3)
    one = F.run_scan(symbols, coordinates, bonds, quad, step=60.0, relaxed=True, runner=runner)
    two = F.run_scan(symbols, coordinates, bonds, quad, step=60.0, relaxed=True, both_directions=True, runner=runner)
    assert len(two["e_qm"]) == 6
    # . two independent xtb runs (multithreaded: not bit-reproducible) -> small tolerance
    assert all(b <= a + 0.05 for a, b in zip(one["e_qm"], two["e_qm"]))
    # . the relaxed barrier is lower than the rigid one
    assert max(one["e_qm"]) - min(one["e_qm"]) < max(rigid["e_qm"]) - min(rigid["e_qm"])
