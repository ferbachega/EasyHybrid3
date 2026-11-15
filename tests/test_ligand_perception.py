"""
Percepcao quimica nativa de ligantes (src/pdynamo/ligand_perception.py), que
substitui o "obabel -p pH" na preparacao de sistemas (2026-10-14).

  1) a partir SO dos atomos pesados (coordenadas 3D geradas pelo RDKit, como
     um ligante cristalografico sem hidrogenios), a molecula reconstruida a
     pH 7.4 -- ordens de ligacao, hidrogenios e cargas formais -- e' a esperada;
  2) com hidrogenios na entrada eles sao mantidos e nada e' protonado;
  3) opls/ligand.ligand_chemistry usa o caminho nativo por padrao, sem obabel.

A validacao completa (533 ligantes de referencia + 97 ligantes cristalograficos
do PLINDER contra o CCD) esta' em DYFF_validation/scripts/ligprep_*.py.
"""
import os
import sys

import numpy as np
import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_TESTS_DIR)
for p in (os.path.join(_REPO, "src"), os.path.join(_REPO, "src", "graphics_engine", "src")):
    if p not in sys.path:
        sys.path.insert(0, p)
pytest.importorskip("rdkit")

from rdkit import Chem, RDLogger                    # noqa: E402
from rdkit.Chem import AllChem                      # noqa: E402
from pdynamo import ligand_perception as LP         # noqa: E402

RDLogger.DisableLog("rdApp.*")
BT = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}

# . (input molecule, expected form at pH 7.4)
CASES = [
    ("OC(=O)c1ccccc1",                 "O=C([O-])c1ccccc1"),          # carboxylic acid
    ("CCN",                            "CC[NH3+]"),                   # primary amine
    ("CN(C)CCc1ccccc1",                "C[NH+](C)CCc1ccccc1"),        # tertiary amine
    ("NC(=N)c1ccccc1",                 "NC(=[NH2+])c1ccccc1"),        # amidine
    ("c1ccncc1",                       "c1ccncc1"),                   # pyridine: neutral
    ("CC(=O)Nc1ccccc1",                "CC(=O)Nc1ccccc1"),            # amide
    ("Oc1ccccc1",                      "Oc1ccccc1"),                  # phenol
    ("NS(=O)(=O)c1ccccc1",             "NS(=O)(=O)c1ccccc1"),         # sulfonamide: neutral
    ("Cn1cnc2c(N)ncnc21",              "Cn1cnc2c(N)ncnc21"),          # 9-methyladenine (amino form)
    ("Cn1ccc(=O)[nH]c1=O",             "Cn1ccc(=O)[nH]c1=O"),         # 1-methyluracil (lactam)
    ("CC(=O)c1ccc(Cl)cc1",             "CC(=O)c1ccc(Cl)cc1"),         # ketone + halogen
    ("Cc1ccc(C#N)s1",                  "Cc1ccc(C#N)s1"),              # thiophene + nitrile
    ("OCC1OC(O)C(O)C(O)C1O",           "OCC1OC(O)C(O)C(O)C1O"),       # glucose
    ("O=[N+]([O-])c1ccccc1",           "O=[N+]([O-])c1ccccc1"),       # nitro keeps its separation
    ("COP(=O)(O)O",                    "COP(=O)([O-])[O-]"),          # phosphate monoester
]


def _heavy_geometry(smiles, seed=7):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    return mol, [mol.GetAtomWithIdx(i).GetSymbol() for i in heavy], [list(conf.GetAtomPosition(i)) for i in heavy]


def _smiles(symbols, charges, bonds):
    em = Chem.RWMol()
    for s, q in zip(symbols, charges):
        a = Chem.Atom(s.capitalize()); a.SetFormalCharge(int(q)); a.SetNoImplicit(True); em.AddAtom(a)
    for i, j, o in bonds:
        em.AddBond(int(i), int(j), BT[int(o)])
    m = em.GetMol()
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=False)


def _canonical(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=False)


@pytest.mark.parametrize("smiles, expected", CASES)
def test_heavy_atoms_to_ph74_form(smiles, expected):
    _, sym, xyz = _heavy_geometry(smiles)
    r = LP.perceive_ligand(sym, xyz, ph=7.4)
    assert _smiles(r["symbols"], r["formal_charges"], r["bonds"]) == _canonical(expected)
    assert r["total_charge"] == Chem.GetFormalCharge(Chem.MolFromSmiles(expected))
    # . every added hydrogen sits ~1 A from its heavy atom
    for heavy, h in r["added_hydrogens"]:
        assert 0.9 < np.linalg.norm(r["coords"][heavy] - r["coords"][h]) < 1.15


def test_hydrogens_in_the_input_are_kept():
    mol, _, _ = _heavy_geometry("CCN")
    conf = mol.GetConformer()
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    r = LP.perceive_ligand(sym, xyz, ph=7.4)
    assert r["added_hydrogens"] == [] and r["decisions"] == [] and r["total_charge"] == 0
    assert len(r["symbols"]) == mol.GetNumAtoms()


def test_ph_controls_protonation():
    _, sym, xyz = _heavy_geometry("OC(=O)c1ccccc1")
    assert LP.perceive_ligand(sym, xyz, ph=2.0)["total_charge"] == 0
    assert LP.perceive_ligand(sym, xyz, ph=7.4)["total_charge"] == -1


def test_ligand_chemistry_is_native_by_default(monkeypatch):
    from pdynamo.opls import ligand as L
    monkeypatch.delenv("EASYHYBRID_LIGAND_CHEMISTRY", raising=False)
    monkeypatch.setenv("PATH", "")                  # . no obabel reachable
    _, sym, xyz = _heavy_geometry("NC(=N)c1ccccc1")
    atoms = [("{}{}".format(s, k + 1), L.ELEMENT_Z[s], tuple(p)) for k, (s, p) in enumerate(zip(sym, xyz))]
    lig = L.ligand_chemistry("BAM", atoms, ph=7.4)
    assert lig.total_charge == 1
    assert sum(1 for (_, z, _) in lig.atoms if z == 1) == 9
