"""
Biblioteca de torcoes DFT por ambiente quimico (src/pdynamo/bond_library.py).

  1) os codigos de ambiente (r1, r2, exato) sao IDENTICOS aos do fraglib
     (fraglib/fraglib/environment.py), com a molecula montada como o EasyHybrid
     a monta (numeros atomicos + ordens Kekule + cargas formais + 3D);
  2) busca exata -> r2 -> r1;
  3) refine_torsions com uma biblioteca cobrindo todas as ligacoes nao roda
     nenhuma varredura xTB e guarda as series da biblioteca.
"""
import os
import sys

import numpy as np
import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_TESTS_DIR)
for p in (os.path.join(_REPO, "src"), os.path.join(_REPO, "fraglib")):
    if p not in sys.path:
        sys.path.insert(0, p)
pytest.importorskip("rdkit")

from pdynamo import bond_library as BL        # noqa: E402

SMILES = ["CC(=O)Nc1ccccc1", "OC[C@H](N)C(=O)[O-]", "C/C=C/C(=O)OCC", "c1ccc(-c2ccncc2)cc1", "C[N+](C)(C)CCO", "CS(=O)(=O)Nc1ccc(Cl)cc1"]


def _molecule(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=3)
    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule, clearAromaticFlags=False)
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondTypeAsDouble()), b.GetIsAromatic()) for b in kekule.GetBonds()]
    return mol, bonds, mol.GetConformer().GetPositions()


@pytest.mark.parametrize("smiles", SMILES)
def test_codes_identical_to_fraglib(smiles):
    from rdkit import Chem
    from fraglib.environment import term_codes, mol_from_mapped_smiles
    mol, bonds, xyz = _molecule(smiles)
    mine = BL.rdkit_molecule([a.GetAtomicNum() for a in mol.GetAtoms()], bonds, [a.GetFormalCharge() for a in mol.GetAtoms()], xyz)
    # . reference built as fraglib builds the library's molecules: mapped SMILES
    #   with explicit hydrogens (AddHs hydrogens are not sanitized: hybridization
    #   UNSPECIFIED instead of S, which changes the codes)
    stereo = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(stereo)
    for a in stereo.GetAtoms():
        a.SetAtomMapNum(a.GetIdx() + 1)
    reference = mol_from_mapped_smiles(Chem.MolToSmiles(stereo, allHsExplicit=True))
    smi = Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True)
    for a, b, order, aromatic in bonds:
        ref = term_codes(reference, (a, b), "bond", smiles=smi)
        got = BL.bond_codes(mine, a, b)
        assert (got["exact"], got["r1"], got["r2"]) == (ref["exact"], ref["env_r1"], ref["env_r2"])


def test_lookup_order():
    mol, bonds, xyz = _molecule("CC(=O)Nc1ccccc1")
    z = [a.GetAtomicNum() for a in mol.GetAtoms()]
    rd = BL.rdkit_molecule(z, bonds, None, xyz)
    j, k = 3, 4                                                    # N-c(aryl)
    codes = BL.bond_codes(rd, j, k)
    library = {"levels": {"r1": {codes["r1"]: {"pairs": [[1.0, 2]], "n": 7}},
                          "r2": {codes["r2"]: {"pairs": [[2.0, 2]], "n": 3}}}}
    assert BL.library_torsions(z, bonds, None, [(j, k)], xyz, library)[(j, k)] == ("r2", [(2.0, 2)], 3)
    del library["levels"]["r2"]
    assert BL.library_torsions(z, bonds, None, [(j, k)], xyz, library)[(j, k)] == ("r1", [(1.0, 2)], 7)
    assert BL.library_torsions(z, bonds, None, [(j, k)], xyz, library, levels=("exact", "r2")) == {}


def test_refine_torsions_uses_library_without_scans(make_vismol_object, monkeypatch):
    pytest.importorskip("pBabel")
    from gui.windows.builder import dyff_parameters as D
    from pdynamo import torsion_fit as TF
    mol, bonds, xyz = _molecule("CCCO")
    v = make_vismol_object([a.GetSymbol() for a in mol.GetAtoms()], [x for b in bonds for x in b[:2]])
    for k, p in enumerate(xyz):
        v.frames[0, k] = p
    v.manual_bonds = {(min(a, b), max(a, b)) for a, b, o, ar in bonds}
    v.manual_bond_orders = {b: 1 for b in v.manual_bonds}
    v.manual_aromatic_bonds = set()
    targets = D.rotatable_bonds(v)
    rd = BL.rdkit_molecule([a.GetAtomicNum() for a in mol.GetAtoms()], bonds, None, xyz)
    library = {"env_code_version": BL.ENV_CODE_VERSION, "levels": {"r2": {
        BL.bond_codes(rd, *bond)["r2"]: {"pairs": [[0.5 + i, 3]], "n": 2} for i, (bond, quad) in enumerate(targets)}}}
    monkeypatch.setattr(BL, "load_library", lambda path=None, charged=False: library)
    monkeypatch.setattr(TF, "run_scan", lambda *a, **k: pytest.fail("no xTB scan expected"))
    results = D.refine_torsions(v, 0, 1, None, workers=1)
    assert len(results) == len(targets) and all(r["source"] == "library" and r["level"] == "r2" for r in results)
    assert D.store_torsion_edits(v, results) == len(targets)
    assert v.dyff_parameter_edits["dihedral_bonds"]["{}-{}".format(*targets[0][0])] == [[0.5, 3]]


@pytest.mark.parametrize("smiles", SMILES)
def test_angle_codes_identical_to_fraglib(smiles):
    from rdkit import Chem
    from fraglib.environment import term_codes, mol_from_mapped_smiles
    mol, bonds, xyz = _molecule(smiles)
    mine = BL.rdkit_molecule([a.GetAtomicNum() for a in mol.GetAtoms()], bonds, [a.GetFormalCharge() for a in mol.GetAtoms()], xyz)
    stereo = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(stereo)
    for a in stereo.GetAtoms():
        a.SetAtomMapNum(a.GetIdx() + 1)
    reference = mol_from_mapped_smiles(Chem.MolToSmiles(stereo, allHsExplicit=True))
    smi = Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True)
    for b in mol.GetAtoms():
        nb = [n.GetIdx() for n in b.GetNeighbors()]
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                t = (min(nb[i], nb[j]), b.GetIdx(), max(nb[i], nb[j]))
                ref = term_codes(reference, t, "angle", smiles=smi)
                got = BL.term_codes(mine, t, "angle")
                assert (got["exact"], got["r1"], got["r2"]) == (ref["exact"], ref["env_r1"], ref["env_r2"])


def test_hessian_calibration_applies_shift_and_ratio():
    from pdynamo import seminario as S
    mol, bonds, xyz = _molecule("CCO")
    z = [a.GetAtomicNum() for a in mol.GetAtoms()]
    rd = BL.rdkit_molecule(z, bonds, None, xyz)
    code = BL.term_codes(rd, (0, 1), "bond")["r2"]
    table = {"table": {"bond": {"r2": {code: [0.01, 2.0, 5]}}, "angle": {}}}
    cal = BL.HessianCalibration(z, bonds, None, xyz, table)
    assert cal("bond", (0, 1)) == (0.01, 2.0, "r2") and cal("bond", (1, 2)) is None
    pairs = [(a, b) for a, b, o, ar in bonds]
    h = np.zeros((3 * len(z), 3 * len(z)))                # springs along every bond (Seminario reads -H_ab)
    for a, b in pairs:
        for i, j, sign in ((a, a, 1), (b, b, 1), (a, b, -1), (b, a, -1)):
            h[3 * i:3 * i + 3, 3 * j:3 * j + 3] += sign * 2000.0 * np.eye(3)
    plain = S.dyff_edits_from_hessian(h, xyz, pairs)
    calibrated = S.dyff_edits_from_hessian(h, xyz, pairs, calibration=cal)
    (b0, k), (b1, k1) = plain["bond_terms"]["0-1"], calibrated["bond_terms"]["0-1"]
    assert b1 == pytest.approx(b0 + 0.01) and k1 == pytest.approx(2.0 * k)
    assert calibrated["bond_terms"]["1-2"] == plain["bond_terms"]["1-2"]
