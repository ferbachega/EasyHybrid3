"""
Double/triple bonds in the "lines" representation (2026-10-14).

vismol.libgl.representations.multiple_bond_split separates the single bonds
(first pass, GL_LINES) from the multiple ones (second pass "lines_multi",
GL_LINES_ADJACENCY, 4 indexes per bond: ref, i, j, j). For a double bond in a
ring, ref must be a ring neighbour of i (the inner line is drawn on its side).
"""
import os
import sys

import pytest

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (os.path.join(_REPO, "src", "graphics_engine", "src"),):
    if p not in sys.path:
        sys.path.insert(0, p)

pytest.importorskip("OpenGL")
from vismol.libgl import representations as R      # noqa: E402


def test_benzaldehyde_nitrile():
    # . Kekule benzene ring 0-5, exocyclic C6=O7, H on C0 (8), separate C9#N10
    pairs = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (5, 6), (6, 7), (0, 8), (9, 10)]
    orders = [2, 1, 2, 1, 2, 1, 1, 2, 1, 3]
    singles, groups = R.multiple_bond_split(pairs, orders, pairs, is_hydrogen=lambda i: i == 8)
    assert sorted(map(tuple, singles.reshape(-1, 2).tolist())) == [(0, 8), (1, 2), (3, 4), (5, 0), (5, 6)]
    ring = groups[R.LINES_MULTI_RING_DOUBLE].reshape(-1, 4).tolist()
    assert ring == [[5, 0, 1, 1], [1, 2, 3, 3], [3, 4, 5, 5]]       # . ref = ring neighbour, never the H
    assert groups[R.LINES_MULTI_DOUBLE].reshape(-1, 4).tolist() == [[5, 6, 7, 7]]   # . ref = substituent (plane)
    assert groups[R.LINES_MULTI_TRIPLE].reshape(-1, 4)[:, 1:3].tolist() == [[9, 10]]


def test_only_singles_and_clamping():
    pairs = [(0, 1), (1, 2)]
    singles, groups = R.multiple_bond_split(pairs, [1, 1], pairs)
    assert groups == {} and singles.tolist() == [0, 1, 1, 2]
    # . orders > 3 are drawn as triple, < 1 as single
    singles, groups = R.multiple_bond_split(pairs, [0, 5], pairs)
    assert singles.tolist() == [0, 1] and list(groups) == [R.LINES_MULTI_TRIPLE]


def test_subset_uses_all_bonds_for_rings():
    # . the representation shows only the double bond 0=1, but the ring 0-1-2-3-4-5
    #   is known from all bonds of the object
    ring = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    singles, groups = R.multiple_bond_split([(0, 1)], [2], ring)
    assert singles.size == 0
    assert groups[R.LINES_MULTI_RING_DOUBLE].tolist() == [5, 0, 1, 1]


def test_sticks_side_vectors():
    # . ethylene-like C0=C1 with substituents in the xy plane, and a C2#N3 triple:
    #   the double bond's side vector lies in that plane, perpendicular to the bond
    import numpy as np
    pairs = [(0, 1), (0, 4), (1, 5), (2, 3)]
    coords = np.array([[0, 0, 0], [1.34, 0, 0], [5, 0, 0], [6.16, 0, 0], [-0.6, 0.9, 0], [1.9, 0.9, 0]], dtype=float)
    mode, ref = R.multiple_bond_references(pairs, [2, 1, 1, 3], pairs)
    assert mode.tolist() == [R.LINES_MULTI_DOUBLE, -1, -1, R.LINES_MULTI_TRIPLE] and ref[0] in (4, 5)
    sides = R.multiple_bond_sides(pairs, mode, ref, coords)
    assert np.allclose(sides[0], [0, 1, 0, R.LINES_MULTI_DOUBLE], atol=1e-6)
    assert sides[3, 3] == R.STICKS_TRIPLE_TRIANGLE and abs(sides[3, 0]) < 1e-6   # . triple: triangle, perpendicular to the bond
