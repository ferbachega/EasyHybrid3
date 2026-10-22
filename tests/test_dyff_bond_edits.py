"""
Edicoes de diedro POR LIGACAO no DYFF (src/pdynamo/dyff_edits.py,
parameterEdits["dihedral_bonds"]): usadas pela biblioteca de fragmentos
(fraglib) para aplicar uma correcao ajustada a todos os diedros em volta de
uma ligacao, independente dos tipos.

  E(com a edicao) = E(diedros da ligacao zerados) + sum_d sum_n c_n cos(n phi_d)
  e os demais diedros (outras ligacoes, mesmo tipo) nao mudam.
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
pytest.importorskip("pBabel")
pytest.importorskip("rdkit")
os.environ.setdefault("PDYNAMO3_PARAMETERS", os.path.join(_PDYNAMO, "parameters"))


def _system_and_coordinates(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from fraglib.mm_baseline import pdynamo_system
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=3)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(a.GetIdx() + 1)
    mapped = Chem.MolToSmiles(mol, allHsExplicit=True)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    system, rd = pdynamo_system(mapped, 0)
    return system, mol, mol.GetConformer().GetPositions()


def _energy(system, xyz, edits):
    from pScientific.Geometry3 import Coordinates3
    from pMolecule.NBModel import NBModelFull
    from pdynamo.dyff_edits import dyff_model
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        system.DefineMMModel(dyff_model("dyff-1.0", edits), log=None)
        system.DefineNBModel(NBModelFull.WithDefaults())
        c3 = Coordinates3.WithExtent(len(system.atoms))
        for a, p in enumerate(xyz):
            for k in range(3):
                c3[a, k] = float(p[k])
        system.coordinates3 = c3
        return float(system.Energy(log=None))


def _dihedral(x, q):
    b0, b1, b2 = x[q[0]] - x[q[1]], x[q[2]] - x[q[1]], x[q[3]] - x[q[2]]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w))


def test_bond_edit_replaces_only_the_dihedrals_around_that_bond():
    system, mol, xyz = _system_and_coordinates("CCCCC")          # pentane: two equivalent C-C-C-C types
    # . the bond C1-C2 (atoms 1, 2)
    j, k = 1, 2
    label = "1-2"
    quads = [(i.GetIdx(), j, k, l.GetIdx()) for i in mol.GetAtomWithIdx(j).GetNeighbors() if i.GetIdx() != k
             for l in mol.GetAtomWithIdx(k).GetNeighbors() if l.GetIdx() != j]
    pairs = [(3.0, 1), (-2.0, 2), (1.5, 3)]
    e_zero = _energy(system, xyz, {"dihedral_bonds": {label: []}})
    e_edit = _energy(system, xyz, {"dihedral_bonds": {label: pairs}})
    expected = sum(c * math.cos(n * _dihedral(xyz, q)) for q in quads for c, n in pairs)
    assert abs((e_edit - e_zero) - expected) < 1e-6
    # . the other bonds keep DYFF: zeroing C2-C3 changes the energy differently
    e_plain = _energy(system, xyz, None)
    e_other = _energy(system, xyz, {"dihedral_bonds": {"2-3": []}})
    assert abs(e_plain - e_zero) > 1e-6 and abs(e_plain - e_other) > 1e-6 and abs(e_zero - e_other) > 1e-6


def test_soften_polar_hydrogens():
    """ Only H bonded to N/O/S get the small LJ type (full and 1-4 tables,
        DYFF 1-4 scale); idempotent; LJ energy changes, bonded terms do not. """
    from fraglib import mm_tools as M
    from pdynamo.dyff_edits import soften_polar_hydrogens, polar_hydrogens, POLAR_H_SIGMA, POLAR_H_EPSILON
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("OCC(=O)N"))           # glycolamide: 1 O-H, 2 N-H, 2 C-H
    AllChem.EmbedMolecule(mol, randomSeed=1)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(a.GetIdx() + 1)
    xyz = mol.GetConformer().GetPositions()
    system, _ = M.dyff_system(Chem.MolToSmiles(mol, allHsExplicit=True), 0)
    polar = polar_hydrogens(system)
    assert len(polar) == 3
    before = (M.energy(system, xyz), dict(system.scratch.energyTerms))
    assert soften_polar_hydrogens(system) == 3
    assert soften_polar_hydrogens(system) == 0
    state = system.mmState
    full, one4 = state.ljParameters.__getstate__(), state.ljParameters14.__getstate__()
    t = state.ljTypeIndices[polar[0]]
    assert full["sigmas"][t] == pytest.approx(POLAR_H_SIGMA) and full["epsilons"][t] == pytest.approx(POLAR_H_EPSILON)
    assert one4["epsilons"][t] == pytest.approx(0.5 * POLAR_H_EPSILON)
    nonpolar = [i for i, a in enumerate(system.atoms) if a.atomicNumber == 1 and i not in polar]
    assert all(state.ljTypeIndices[i] != t for i in nonpolar)
    M.energy(system, xyz)
    after = dict(system.scratch.energyTerms)
    assert after["MM/MM Lennard-Jones"] != pytest.approx(before[1]["MM/MM Lennard-Jones"])
    assert after["Cosine Dihedral"] == pytest.approx(before[1]["Cosine Dihedral"])


def test_polar_hydrogens_include_halogen_donors():
    """ H bonded to F/Cl/Br/I is polar (DYFF_validation Phase IV); H on C bonded to a
        halogen-bearing carbon is not. """
    from fraglib import mm_tools as M
    from pdynamo.dyff_edits import polar_hydrogens
    from rdkit import Chem
    for smiles, n in (("Cl", 1), ("Br", 1), ("F", 1), ("I", 1), ("OCCl", 1), ("FC(F)F", 0)):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        for a in mol.GetAtoms():
            a.SetAtomMapNum(a.GetIdx() + 1)
        system, _ = M.dyff_system(Chem.MolToSmiles(mol, allHsExplicit=True), 0)
        assert len(polar_hydrogens(system)) == n, smiles
