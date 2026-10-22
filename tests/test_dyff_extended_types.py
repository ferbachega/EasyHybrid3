"""
EasyHybrid's DYFF parameter set dyff-1.0-eh (src/pdynamo/data/forcefields/dyff/dyff-1.0-eh):
pDynamo3 dyff-1.0 plus ADDED types/patterns.

  1) molecules dyff-1.0 cannot type get the expected types: phosphinine P:Res, borazine B:Tri,
     selenophene / 1,2-selenazole Se:Res, P4 cage P:Tet:+3 (pyramidal, not P:Res), diborane
     bridging H:::b;
  2) nothing that dyff-1.0 types changes type (a sample of common chemistry, atom by atom);
  3) dyff_model("dyff-1.0") uses the extension, EASYHYBRID_DYFF_ORIGINAL=1 restores pDynamo's set.
"""
import contextlib
import io
import os
import sys

import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_TESTS_DIR)
for p in (os.path.join(_REPO, "src"), os.path.join(_REPO, "src", "graphics_engine", "src")):
    if p not in sys.path:
        sys.path.insert(0, p)
pytest.importorskip("pBabel")
Chem = pytest.importorskip("rdkit.Chem")


def _system(mol, bonds=None):
    import pBabel  # noqa: F401
    from pMolecule import Atom, Bond, BondType, Connectivity, ConvertInputConnectivity, System
    c = Connectivity()
    for a in mol.GetAtoms():
        c.AddNode(Atom.WithOptions(atomicNumber=a.GetAtomicNum(), label="{}{}".format(a.GetSymbol(), a.GetIdx() + 1),
                                   formalCharge=a.GetFormalCharge()))
    types = {1.0: BondType.Single, 2.0: BondType.Double, 3.0: BondType.Triple}
    if bonds is None:
        k = Chem.Mol(mol); Chem.Kekulize(k, clearAromaticFlags=False)
        bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble(), b.GetIsAromatic()) for b in k.GetBonds()]
    for i, j, order, aromatic in bonds:
        c.AddEdge(Bond.WithNodes(c.nodes[i], c.nodes[j], isAromatic=aromatic, type=types.get(order, BondType.Single)))
    ConvertInputConnectivity(c, {})
    return System.FromConnectivity(connectivity=c)


def _types(system, original=False):
    from pdynamo.dyff_edits import dyff_model
    old = os.environ.pop("EASYHYBRID_DYFF_ORIGINAL", None)
    if original: os.environ["EASYHYBRID_DYFF_ORIGINAL"] = "1"
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            system.DefineMMModel(dyff_model("dyff-1.0", None), log=None)
    finally:
        os.environ.pop("EASYHYBRID_DYFF_ORIGINAL", None)
        if old is not None: os.environ["EASYHYBRID_DYFF_ORIGINAL"] = old
    idx = list(system.mmState.atomTypeIndices); names = list(system.mmState.atomTypes)
    return [names[i] for i in idx]


def _smiles(smi):
    return _system(Chem.AddHs(Chem.MolFromSmiles(smi)))


@pytest.mark.parametrize("smi, element, label", [
    ("c1ccpcc1", 15, "P:Res"), ("B1NBNBN1", 5, "B:Tri"), ("c1cc[se]c1", 34, "Se:Res"),
    ("c1cn[se]c1", 34, "Se:Res"), ("O=c1cc[se][nH]1", 34, "Se:Res"), ("P12P3P1P23", 15, "P:Tet:+3")])
def test_new_types(smi, element, label):
    s = _smiles(smi)
    types = _types(s)
    assert {t for a, t in zip(s.atoms, types) if a.atomicNumber == element} == {label}


def test_original_set_fails_on_them():
    from pMolecule.MMModel.MMModelError import MMModelError
    for smi in ("c1ccpcc1", "B1NBNBN1", "c1cc[se]c1", "P12P3P1P23"):
        with pytest.raises(MMModelError):
            _types(_smiles(smi), original=True)


def test_diborane_bridges():
    m = Chem.RWMol()
    for e in ["B", "B", "H", "H", "H", "H", "H", "H"]:
        a = Chem.Atom(e); a.SetNoImplicit(True); m.AddAtom(a)
    bonds = [(0, 2, 1.0, False), (0, 3, 1.0, False), (1, 4, 1.0, False), (1, 5, 1.0, False),
             (0, 6, 1.0, False), (1, 6, 1.0, False), (0, 7, 1.0, False), (1, 7, 1.0, False)]
    types = _types(_system(m.GetMol(), bonds))
    assert types[:2] == ["B:Tet", "B:Tet"] and types[6:] == ["H:::b", "H:::b"]


@pytest.mark.parametrize("smi", ["c1ccccc1", "c1ccncc1", "c1ccsc1", "CC(=O)NC", "[NH3+]CC(=O)[O-]", "NS(=O)(=O)c1ccccc1",
                                 "OP(=O)(O)O", "c1ccc2[nH]ccc2c1", "CSC", "O=C1NC(=O)C=C1", "ClC(Cl)Cl", "CC#N", "C[Se]C", "PP"])
def test_unchanged_for_dyff_1_0_molecules(smi):
    assert _types(_smiles(smi)) == _types(_smiles(smi), original=True)


def test_dyff_model_uses_extension():
    import pBabel  # noqa: F401
    from pdynamo.dyff_edits import dyff_model
    assert dyff_model("dyff-1.0").parameterSet == "dyff-1.0-eh"


def _opt_angle(smi, sym, original=False):
    import numpy as np
    from rdkit.Chem import AllChem
    import pBabel  # noqa: F401
    from pMolecule.NBModel import NBModelFull
    from pScientific.Geometry3 import Coordinates3
    from pSimulation import LBFGSMinimize_SystemGeometry
    from pdynamo.dyff_edits import dyff_model
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=3)
    s = _system(m); x = m.GetConformer().GetPositions()
    c3 = Coordinates3.WithExtent(len(x))
    for i in range(len(x)):
        for k in range(3): c3[i, k] = float(x[i][k])
    s.coordinates3 = c3
    old = os.environ.pop("EASYHYBRID_DYFF_ORIGINAL", None)
    if original: os.environ["EASYHYBRID_DYFF_ORIGINAL"] = "1"
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            s.DefineMMModel(dyff_model("dyff-1.0", None), log=None); s.DefineNBModel(NBModelFull.WithDefaults())
            LBFGSMinimize_SystemGeometry(s, log=None, maximumIterations=5000, rmsGradientTolerance=0.001)
    finally:
        os.environ.pop("EASYHYBRID_DYFF_ORIGINAL", None)
        if old is not None: os.environ["EASYHYBRID_DYFF_ORIGINAL"] = old
    y = np.array([[s.coordinates3[i, k] for k in range(3)] for i in range(len(x))])
    a = next(a for a in m.GetAtoms() if a.GetSymbol() == sym); j = a.GetIdx(); i, k = [n.GetIdx() for n in a.GetNeighbors() if n.GetSymbol() != "H"][:2]
    u, v = y[i] - y[j], y[k] - y[j]
    return float(np.degrees(np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v))))


def test_bent_resonant_angles():
    """ P:Res / Se:Res angles use the general form with calibrated natural angles: the parent
        rings reproduce their references; thiophene (S:Res, a dyff-1.0 type) is unchanged. """
    assert abs(_opt_angle("c1cc[se]c1", "Se") - 87.77) < 0.3
    assert abs(_opt_angle("c1ccpcc1", "P") - 99.9) < 0.3
    assert abs(_opt_angle("c1ccsc1", "S") - _opt_angle("c1ccsc1", "S", original=True)) < 1e-6
