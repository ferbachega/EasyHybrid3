"""
Codigos de ambiente sem RDKit (src/pdynamo/env_codes.py), 2026-10-14.

Os invariantes de atomo (aromaticidade, anel, hibridizacao, cargas apos o
"cleanUp" de nitro/P/halogenios) e os codigos r1/r2 de cada ligacao devem ser
IDENTICOS aos do RDKit (as chaves das bibliotecas DFT foram feitas com ele).
A comparacao completa (18 413 moleculas, 0 diferencas) esta' em
DYFF_validation/scripts/env_codes_check.py.
"""
import os
import sys

import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_TESTS_DIR)
for p in (os.path.join(_REPO, "src"),):
    if p not in sys.path:
        sys.path.insert(0, p)

from pdynamo import env_codes as E          # noqa: E402

# . aromatic heterocycles, fused/azulene-like systems, pyridone, nitro (cleanUp),
#   N-oxide, sulfone, phosphonate, perchlorate-like, azide, macrocyclic ether,
#   anions/cations, radical
SMILES = ["c1ccccc1", "c1ccncc1", "c1cc[nH]c1", "c1ccoc1", "c1ccsc1", "O=c1cccc[nH]1", "c1ccc2ccccc2c1",
          "c1ccc2cccc2cc1", "C[N+](=O)[O-]", "CN(=O)=O", "[O-][n+]1ccccc1", "CS(=O)(=O)c1ccccc1",
          "COP(=O)(OC)C", "OCl(=O)(=O)=O", "CN=[N+]=[N-]", "C1COCCOCCOCCOCCO1", "CC(=O)[O-]",
          "C[NH3+]", "Cc1n[nH]c2ccccc12", "O=C1NC(=O)c2ccccc21", "c1ccc2[nH]ccc2c1", "C=CC=O",
          "N#Cc1ccccc1", "[O]P(=O)(O)OC", "Cc1ccc2c(c1)Sc1ccccc1[S@@]2(C)c1ccccc1", "c1cc[se]c1",
          "B1c2ccccc2-c2ccccc12", "C1=CC=CC=CC=C1", "c1cc2ccc3cccc4ccc(c1)c2c34"]


@pytest.mark.parametrize("smiles", SMILES)
def test_identical_to_rdkit(smiles):
    Chem = pytest.importorskip("rdkit.Chem")
    from pdynamo import bond_library as BL
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Chem.SanitizeMol(mol)                    # . as the library / lookup molecules (AddHs leaves the new H unsanitized)
    kek = Chem.Mol(mol)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    bt = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}
    own = E.Molecule([a.GetAtomicNum() for a in kek.GetAtoms()],
                     [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), bt[b.GetBondType()]) for b in kek.GetBonds()],
                     [a.GetFormalCharge() for a in kek.GetAtoms()])
    for a in mol.GetAtoms():
        assert own.atom_invariant(a.GetIdx()) == BL._atom_invariant(a), (smiles, a.GetIdx())
    ids_r = {k: BL.atom_identifiers(mol, k) for k in (1, 2)}
    ids_o = {k: E.atom_identifiers(own, k) for k in (1, 2)}
    for b in mol.GetBonds():
        t = (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        ref = BL.term_codes(mol, t, "bond", identifiers=ids_r, smiles="-", ranks=list(range(mol.GetNumAtoms())))
        assert E.term_codes_r(own, t, "bond", ids_o) == {"r1": ref["r1"], "r2": ref["r2"]}


def test_no_rdkit_needed():
    # . benzene from atomic numbers and Kekule bonds alone
    z = [6] * 6 + [1] * 6
    bonds = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)] + [(i, i + 6, 1) for i in range(6)]
    m = E.Molecule(z, bonds)
    assert m.aromatic_atoms == set(range(6))
    assert all(m.hybridization[i] == "SP2" for i in range(6))
    assert {m.bond_symbol(b) for b in range(6)} == {"a"}


def test_exact_code_is_numbering_invariant():
    # . ethanol in two atom orders: same exact code for the C-O bond
    m1 = E.Molecule([6, 6, 8, 1, 1, 1, 1, 1, 1], [(0, 1, 1), (1, 2, 1), (0, 3, 1), (0, 4, 1), (0, 5, 1), (1, 6, 1), (1, 7, 1), (2, 8, 1)])
    m2 = E.Molecule([8, 6, 6, 1, 1, 1, 1, 1, 1], [(0, 1, 1), (1, 2, 1), (2, 3, 1), (2, 4, 1), (2, 5, 1), (1, 6, 1), (1, 7, 1), (0, 8, 1)])
    assert E.exact_code(m1, (1, 2)) == E.exact_code(m2, (1, 0))
    assert E.exact_code(m1, (0, 1)) != E.exact_code(m1, (1, 2))


def test_cleanup_like_rdkit():
    # . neutral 5-valent forms as a user structure may have them: CN(=O)=O,
    #   CN=N#N, C=P(=O)C, OCl(=O)(=O)=O -> RDKit's charge-separated forms
    cases = {
        "nitro":       ([6, 7, 8, 8], [(0, 1, 1), (1, 2, 2), (1, 3, 2)], {1: 1, 2: -1, 3: 0}),
        "azide":       ([6, 7, 7, 7], [(0, 1, 1), (1, 2, 2), (2, 3, 3)], {2: 1, 3: -1}),
        "ylide":       ([6, 15, 8, 6, 1, 1, 1, 1, 1], [(0, 1, 2), (1, 2, 2), (1, 3, 1), (0, 4, 1), (0, 5, 1), (3, 6, 1), (3, 7, 1), (3, 8, 1)],
                        {1: 1, 2: -1}),                 # . explicit H: the CH2 has degree 3 (RDKit gives C=[P+]([O-])C)
        "perchlorate": ([8, 17, 8, 8, 8], [(0, 1, 1), (1, 2, 2), (1, 3, 2), (1, 4, 2)], {1: 3, 2: -1, 3: -1, 4: -1}),
    }
    for name, (z, bonds, expected) in cases.items():
        m = E.Molecule(z, bonds)
        for a, q in expected.items():
            assert m.charge[a] == q, (name, a, m.charge)


def test_library_lookups_without_rdkit(monkeypatch):
    # . the DFT torsion library and the Hessian calibration with RDKit unimportable
    import importlib
    import builtins
    real_import = builtins.__import__
    def no_rdkit(name, *args, **kwargs):
        if name == "rdkit" or name.startswith("rdkit."):
            raise ImportError("rdkit blocked in this test")
        return real_import(name, *args, **kwargs)
    monkeypatch.delenv("EASYHYBRID_ENV_CODES", raising=False)
    monkeypatch.setattr(builtins, "__import__", no_rdkit)
    from pdynamo import bond_library as BL
    BL = importlib.reload(BL)
    from pdynamo import residue_library as RL
    if BL.load_library() is None or BL.load_hessian_library() is None:
        pytest.skip("library files not installed")
    assert BL.load_library()["env_code_version"] == E.ENV_CODE_VERSION_OWN
    # . N-butylbenzamide, Kekule bonds, all H explicit (heavy: c0-c5 ring, C6(=O7)N8 C9 C10 C11 C12)
    heavy = [6] * 6 + [6, 8, 7, 6, 6, 6, 6]
    bonds = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1), (6, 7, 2), (6, 8, 1),
             (8, 9, 1), (9, 10, 1), (10, 11, 1), (11, 12, 1)]
    nh = {1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 8: 1, 9: 2, 10: 2, 11: 2, 12: 3}
    z = list(heavy)
    for a, n in nh.items():
        for _ in range(n):
            z.append(1); bonds.append((a, len(z) - 1, 1))
    bonds4 = [(i, j, o, False) for i, j, o in bonds]
    xyz = [[0.0, 0.0, 0.0]] * len(z)                    # . coordinates are not used by the codes
    edits, report = RL.library_dihedrals(z, bonds4, [0] * len(z), xyz, list(range(len(z))), charged=False)
    assert report["rotatable"] >= 4 and report["torsions"] >= 3, report
    cal = BL.HessianCalibration(z, bonds4, [0] * len(z), xyz)
    assert cal("bond", (6, 7)) is not None
    nb = {}
    for i, j, _ in bonds:
        nb.setdefault(i, []).append(j); nb.setdefault(j, []).append(i)
    angles = [(a, b, c) for b, ns in nb.items() for k, a in enumerate(ns) for c in ns[k + 1:]]
    assert sum(cal("angle", t) is not None for t in angles) > 0
    monkeypatch.setattr(builtins, "__import__", real_import)
    importlib.reload(BL)
