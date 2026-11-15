"""
Metodo de Seminario modificado (src/pdynamo/seminario.py) e edicoes por
termo de ligacao/angulo do DYFF (dyff_edits "bond_terms"/"angle_terms").

  1) unidades: frequencias recalculadas da Hessiana do xtb == as do xtb;
  2) ida e volta: Hessiana do PROPRIO DYFF (diferencas finitas) no minimo DYFF
     -> Seminario -> constantes de ligacao proximas das do DYFF;
  3) a edicao por termo muda so aquele termo: E(edicao) - E(original) ==
     diferenca analitica do termo.
"""
import math
import os
import sys

import numpy as np
import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_TESTS_DIR)
_PDYNAMO = os.environ.get("PDYNAMO3_HOME", os.path.expanduser("~/programs/pDynamo3"))
for p in (_PDYNAMO, os.path.join(_REPO, "src"), os.path.join(_REPO, "fraglib")):
    if p not in sys.path:
        sys.path.insert(0, p)
os.environ.setdefault("PDYNAMO3_PARAMETERS", os.path.join(_PDYNAMO, "parameters"))

from pdynamo import seminario as S        # noqa: E402


def _mapped(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    xyz = mol.GetConformer().GetPositions()
    for a in mol.GetAtoms():
        a.SetAtomMapNum(a.GetIdx() + 1)
    return Chem.MolToSmiles(mol, allHsExplicit=True), mol, xyz


def test_frequencies_from_xtb_hessian():
    from pdynamo.opls.ligand import find_xtb
    if find_xtb() is None:
        pytest.skip("xtb not found")
    sym = ["O", "H", "H"]
    xyz, h, fx = S.xtb_hessian(sym, [[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]])
    assert np.allclose(S.frequencies(h, sym), fx[-3:], atol=0.5)
    (r0, k), = S.bond_constants(h, xyz, [(0, 1)]).values()
    assert 4000 < k < 5000 and abs(r0 - 0.96) < 0.01


def test_seminario_recovers_dyff_bond_constants():
    pytest.importorskip("pBabel")
    from fraglib import mm_tools as M
    mapped, mol, xyz = _mapped("CCO")                      # ethanol
    system, rd = M.dyff_system(mapped, 0)
    x = M.optimize(system, xyz)
    h = M.hessian(system, x)
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in rd.GetBonds()]
    sem = S.bond_constants(h, x, bonds)
    # . DYFF's own bond parameters (E = k (b - b0)^2 -> 1/2 k' with k' = 2 k)
    state = [c for c in system.mmState.mmTerms if type(c).__name__ == "HarmonicBondContainer"][0].__getstate__()
    dyff = {}
    for (i, j), p in [((t[0], t[1]), t[2]) for t in state["terms"]]:
        dyff[(min(i, j), max(i, j))] = 2.0 * state["parameters"][p][1]
    ratios = [sem[b][1] / dyff[b] for b in sem]
    assert 0.85 < np.median(ratios) < 1.15, ratios


def test_per_term_edits_change_only_that_term():
    pytest.importorskip("pBabel")
    from fraglib import mm_tools as M
    mapped, mol, xyz = _mapped("CCO")
    system, rd = M.dyff_system(mapped, 0)
    e0 = M.energy(system, xyz)
    a, b = sorted((rd.GetBonds()[0].GetBeginAtomIdx(), rd.GetBonds()[0].GetEndAtomIdx()))
    r = float(np.linalg.norm(xyz[a] - xyz[b]))
    state = [c for c in system.mmState.mmTerms if type(c).__name__ == "HarmonicBondContainer"][0].__getstate__()
    old = next(state["parameters"][t[2]] for t in state["terms"] if {t[0], t[1]} == {a, b})
    new = (1.30, 1500.0)
    edited, _ = M.dyff_system(mapped, 0, {"bond_terms": {"{}-{}".format(a, b): list(new)}})
    e1 = M.energy(edited, xyz)
    expected = new[1] * (r - new[0]) ** 2 - old[1] * (r - old[0]) ** 2
    assert abs((e1 - e0) - expected) < 1e-6
    # . an angle term: c-n series of seminario.dyff_angle_pairs, curvature check
    pairs = S.dyff_angle_pairs(math.radians(110.0), 500.0)
    t = np.radians([109.0, 110.0, 111.0])
    e = [sum(c * math.cos(n * x) for c, n in pairs) for x in t]
    curvature = (e[0] - 2 * e[1] + e[2]) / (t[1] - t[0]) ** 2
    assert abs(curvature - 500.0) < 1.0 and e[1] < e[0] and e[1] < e[2]


def test_refine_torsions_with_xtb(make_vismol_object):
    """ Hessian bonds/angles + per-bond torsion fit (xTB scans of the molecule)
        on 1-propanol: every rotatable bond fitted, the self-consistent fit does
        not make the profile worse than DYFF, and the edits are stored with the
        topology signature (dropped when the molecule changes). """
    pytest.importorskip("pBabel")
    from pdynamo.opls.ligand import find_xtb
    if find_xtb() is None:
        pytest.skip("xtb not found")
    import pBabel  # noqa: F401
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from gui.windows.builder import dyff_parameters as D
    from pdynamo.dyff_edits import applicable_edits
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCO"))
    AllChem.EmbedMolecule(mol, randomSeed=5)
    xyz = mol.GetConformer().GetPositions()
    pairs = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    v = make_vismol_object([a.GetSymbol() for a in mol.GetAtoms()], [x for p in pairs for x in p])
    for k, p in enumerate(xyz):
        v.frames[0, k] = p
    v.manual_bonds = {(min(a, b), max(a, b)) for a, b in pairs}
    v.manual_bond_orders = {b: 1 for b in v.manual_bonds}
    v.manual_aromatic_bonds = set()
    bonds = D.rotatable_bonds(v)
    assert [b for b, q in bonds] == [(0, 1), (1, 2), (2, 3)]          # C-C, C-C, C-O
    h = D.compute_hessian_edits(v, 0, 1, threads=2)
    base = {"bond_terms": h["edits"]["bond_terms"], "angle_terms": h["edits"]["angle_terms"]}
    results = D.refine_torsions(v, 0, 1, None, step=30.0, workers=3, base_edits=base, use_library=False)
    assert results and all("pairs" in r for r in results)
    assert all(r["fit_rmse"] <= r["before_rmse"] + 0.3 for r in results)
    n = D.store_torsion_edits(v, results, h)
    assert n == len(results) and v.dyff_parameter_edits.get("term_signature")
    edits, dropped = applicable_edits(v.dyff_parameter_edits, len(v.atoms), v.manual_bonds)
    assert not dropped and edits["dihedral_bonds"]
    edits, dropped = applicable_edits(v.dyff_parameter_edits, len(v.atoms) + 1, v.manual_bonds)
    assert dropped and "dihedral_bonds" not in edits and "bond_terms" not in edits


def test_subsystem_rules():
    """ Capped subsystem around a selection (src/pdynamo/subsystem.py): rings
        enter whole, multiple bonds are never cut, a metal bonded to the region
        is pulled in with its ligands, cuts are capped with H along the bond,
        and per-term edits map back only for terms fully inside the selection. """
    from pdynamo import subsystem as SS
    # . toluene-like ring with a chain: C0-C1-C2-C3(=O4), C3 bonded to ring C5..C10
    z = [6, 6, 6, 6, 8] + [6] * 6
    xyz = np.array([[i * 1.5, 0.0, 0.0] for i in range(5)] + [[6 + np.cos(t), np.sin(t), 0] for t in np.linspace(0, 5 * np.pi / 3, 6)])
    bonds = [(0, 1, 1, False), (1, 2, 1, False), (2, 3, 1, False), (3, 4, 2, False), (3, 5, 1, False)] + \
            [(5 + i, 5 + (i + 1) % 6, 1, True) for i in range(6)]
    sub = SS.build_subsystem(z, xyz, bonds, selection=[0], shell=1)
    assert sorted(sub["atoms"]) == [0, 1] and sub["n_caps"] == 1            # chain cut at C1-C2, one cap
    cap = sub["coordinates"][-1]
    assert abs(np.linalg.norm(cap - xyz[1]) - 1.09) < 1e-6                 # along C1->C2, C-H length
    sub = SS.build_subsystem(z, xyz, bonds, selection=[2], shell=1)
    assert {3, 4} <= set(sub["atoms"])                                      # C=O never cut
    assert not set(range(5, 11)) & set(sub["atoms"])                        # ring not reached: cut at the single C3-C5
    sub = SS.build_subsystem(z, xyz, bonds, selection=[2], shell=2)
    assert set(range(5, 11)) <= set(sub["atoms"])                           # one ring atom reached -> the whole ring
    # . metal: Zn bonded to N; selecting the N pulls the Zn and its other ligand in
    z2 = [7, 30, 8, 6]
    xyz2 = np.array([[0, 0, 0], [2.0, 0, 0], [4.0, 0, 0], [-1.4, 0, 0]], float)
    sub2 = SS.build_subsystem(z2, xyz2, [(0, 1, 1, False), (1, 2, 1, False), (0, 3, 1, False)], selection=[0], shell=0)
    assert {0, 1, 2} <= set(sub2["atoms"])
    # . edits back to the parent: only terms with every atom selected
    edits = {"bond_terms": {"0-1": [1.5, 300.0]}, "angle_terms": {"0-1-2": [[1.0, 0]]}}
    back = SS.map_term_edits_to_parent(edits, [10, 11, 12], selection={10, 11})
    assert list(back["bond_terms"]) == ["10-11"] and not back["angle_terms"]


def test_outofplane_round_trip():
    """ "QM" = DYFF's own Hessian with every out-of-plane series doubled: the
        curvature fit recovers scale 2 per centre; the per-centre edit scales
        exactly the out-of-plane energy. """
    pytest.importorskip("pBabel")
    from fraglib import mm_tools as M
    mapped, mol, xyz = _mapped("CC(=O)Nc1ccccc1")          # acetanilide
    make = lambda edits: M.dyff_system(mapped, 0, edits)[0]
    system = make(None)
    x = M.optimize(system, xyz)
    centres = S.oop_centres(system)
    assert len(centres) == 8                                # C=O carbon, N, six ring carbons
    doubled = {"outofplane_centres": {str(j): [[2 * c, n] for c, n in s] for j, (_, s) in centres.items()}}
    target = make(doubled)
    # . the edit doubles the out-of-plane energy and nothing else
    x_bent = x.copy(); x_bent[list(centres)[0]] += 0.15
    terms = []
    for s in (system, target):
        M.energy(s, x_bent); terms.append(dict(s.scratch.energyTerms))
    assert terms[1]["Cosine Out-Of-Plane"] == pytest.approx(2 * terms[0]["Cosine Out-Of-Plane"], rel=1e-9)
    assert terms[1]["Cosine Dihedral"] == pytest.approx(terms[0]["Cosine Dihedral"], rel=1e-9)
    report = {}
    edits = S.outofplane_edits(S.mm_hessian(target, x), x, make, None, report=report)
    for j, r in report.items():
        assert r["scale"] == pytest.approx(2.0, abs=0.05), (j, r)
    assert set(edits["outofplane_centres"]) == {str(j) for j in centres}


def test_subsystem_cuts_macrocycles_not_small_rings():
    """ A 6-ring is never cut; a 16-membered macrocycle (as the disulfide-closed
        backbone of a protein) is cut and capped like a chain. """
    from pdynamo import subsystem as SS
    n = 16 + 6
    bonds = [(i, (i + 1) % 16, 1, False) for i in range(16)]                   # macrocycle 0..15
    bonds += [(16 + i, 16 + (i + 1) % 6, 1, False) for i in range(6)]          # cyclohexane 16..21
    bonds += [(8, 16, 1, False)]
    z = [6] * n
    xyz = np.random.default_rng(1).normal(size=(n, 3)) * 3.0
    assert (0, 15) in SS.ring_bonds(n, bonds) and (0, 15) not in SS.ring_bonds(n, bonds, max_size=8)
    assert (16, 17) in SS.ring_bonds(n, bonds, max_size=8)
    sub = SS.build_subsystem(z, xyz, bonds, [0], shell=1)
    assert sub["atoms"] == [0, 1, 15]                                          # macrocycle cut
    sub = SS.build_subsystem(z, xyz, bonds, [8], shell=1)
    assert set(range(16, 22)) <= set(sub["atoms"])                            # the 6-ring stays whole
