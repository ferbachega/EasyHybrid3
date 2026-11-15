"""
Bond orders of residues without hydrogens (raw PDB), 2026-10-15
(src/graphics_engine/src/vismol/core/residue_bond_orders.py).

Full check on the 97 raw PLINDER complexes:
DYFF_validation/scripts/raw_pdb_bond_orders_check.py.
"""
import os
import sys
from types import SimpleNamespace

import numpy as np

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (os.path.join(_REPO, "src", "graphics_engine", "src"), os.path.join(_REPO, "src")):
    if p not in sys.path:
        sys.path.insert(0, p)

from vismol.core import residue_bond_orders as RBO        # noqa: E402


def _atoms(spec):
    """ spec: [(resn, resi, name, element)] -> {index: atom} """
    residues, atoms = {}, {}
    for k, (resn, resi, name, el) in enumerate(spec):
        res = residues.setdefault((resn, resi), SimpleNamespace(name=resn))
        atoms[k] = SimpleNamespace(symbol=el, name=name, residue=res)
    return atoms


def test_templates():
    assert RBO.template("PHE", {"C", "O"})[("C", "O")] == 2
    assert RBO.template("HOH", {"O"}) == {}
    assert RBO.template("NALA", {"C", "O"}) == {("C", "O"): 2}                 # . Amber N-terminal
    assert RBO.template("HSD", set())[("CE1", "NE2")] == 2                       # . CHARMM HIS(delta)
    assert RBO.template("DA5", set())[("C6", "N1")] == 2                         # . Amber DNA
    assert RBO.template("TEP", set()) is None                                    # . ligand


def test_peptide_without_h():
    # . Gly-Ser heavy atoms: the valence heuristic would see unsaturation everywhere
    spec = [("GLY", 1, "N", "N"), ("GLY", 1, "CA", "C"), ("GLY", 1, "C", "C"), ("GLY", 1, "O", "O"),
            ("SER", 2, "N", "N"), ("SER", 2, "CA", "C"), ("SER", 2, "C", "C"), ("SER", 2, "O", "O"),
            ("SER", 2, "CB", "C"), ("SER", 2, "OG", "O")]
    atoms = _atoms(spec)
    pairs = [0, 1, 1, 2, 2, 3, 2, 4, 4, 5, 5, 6, 6, 7, 5, 8, 8, 9]
    wrong = np.full(9, 2, dtype=np.uint32)                                       # . everything double
    out = RBO.residue_aware_orders(atoms, pairs, wrong, np.zeros((10, 3)))
    assert out.tolist() == [1, 1, 2, 1, 1, 1, 2, 1, 1]                           # . only the two C=O


def test_residues_with_h_are_untouched():
    spec = [("ALA", 1, "C", "C"), ("ALA", 1, "O", "O"), ("ALA", 1, "H", "H")]
    atoms = _atoms(spec)
    orders = np.array([1, 1], dtype=np.uint32)
    assert RBO.residue_aware_orders(atoms, [0, 1, 0, 2], orders, np.zeros((3, 3))).tolist() == [1, 1]


def test_atoms_without_residue_are_untouched():
    atoms = {0: SimpleNamespace(symbol="C", name="C1", residue=None), 1: SimpleNamespace(symbol="O", name="O1", residue=None)}
    assert RBO.residue_aware_orders(atoms, [0, 1], np.array([2], dtype=np.uint32), np.zeros((2, 3))).tolist() == [2]


def test_ligand_by_native_perception():
    # . acetic acid heavy atoms, crystal-like geometry, one residue "ACT"
    xyz = np.array([[0.0, 0.0, 0.0], [1.50, 0.0, 0.0], [2.12, 1.06, 0.0], [2.14, -1.10, 0.0]])
    atoms = _atoms([("ACT", 1, "CH3", "C"), ("ACT", 1, "C", "C"), ("ACT", 1, "O", "O"), ("ACT", 1, "OXT", "O")])
    out = RBO.residue_aware_orders(atoms, [0, 1, 1, 2, 1, 3], np.ones(3, dtype=np.uint32), xyz)
    assert sorted(out.tolist()) == [1, 1, 2] and out[0] == 1                     # . one C=O, C-C single
