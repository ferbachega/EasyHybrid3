"""
Planar conjugated groups (2026-10-15; user report: Gln/Arg side-chain nitrogens of a
prepared protein, 4oas, out of plane).

1. dyff_prep.planarize_conjugated_nh2: the hydrogens the builder puts on an NH2
   bonded to a trigonal carbon go into that carbon's plane (amides, Arg).
2. The force-matched residue library keeps the calibrated torsions of every bond
   between two trigonal atoms (amide, guanidinium, omega): the force matching had
   made some of them non-planar (fraglib/fraglib/fix_planar_torsions.py).
"""
import json
import math
import os
import sys

import pytest

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for p in (os.path.join(_REPO, "src"), os.path.join(_REPO, "fraglib")):
    if p not in sys.path:
        sys.path.insert(0, p)

DATA = os.path.join(_REPO, "src", "pdynamo", "data")


def _even_profile(pairs):
    ev = lambda phi: sum(c * math.cos(n * math.radians(phi)) for c, n in pairs if int(n) % 2 == 0)
    grid = list(range(0, 181, 5))
    e = [ev(a) for a in grid]
    return ev(90) - min(ev(0), ev(180)), grid[e.index(min(e))]


def test_fm_library_planar_bonds_are_planar():
    # . every bond between two trigonal atoms: a single 2-fold term, minimum in the
    #   plane, with the barrier of the calibrated series (fix_planar_torsions.py)
    from fraglib.run_force_matching import is_planar_bond
    fm = json.load(open(os.path.join(DATA, "dyff_residue_library_fm.json")))["templates"]
    cal = json.load(open(os.path.join(DATA, "dyff_residue_library_calibrated.json")))["templates"]
    n = 0
    for name, t in fm.items():
        cmap = {frozenset((a, b)): p for a, b, p in cal[name]["dihedral_bonds"]}
        for a, b, pairs in t["dihedral_bonds"]:
            if is_planar_bond(t, a, b) and frozenset((a, b)) in cmap:
                assert len(pairs) == 1 and int(pairs[0][1]) == 2 and pairs[0][0] < 0, (name, a, b, pairs)
                barrier, minimum = _even_profile(cmap[frozenset((a, b))])
                assert abs(-2.0 * pairs[0][0] - barrier) < 1e-6, (name, a, b)
                n += 1
    assert n == 61
    # . Arg: the three guanidinium C-N bonds have a planar series, CZ-NH2 = CZ-NH1
    arg = fm["ARG"]
    cz = [k for k, at in enumerate(arg["atoms"]) if at["element"] == "C"
          and sum(arg["atoms"][q]["element"] == "N" for q in at["neighbors"]) == 3][0]
    series = {frozenset((a, b)): p for a, b, p in arg["dihedral_bonds"] if "r:%d" % cz in (a, b)}
    assert len(series) == 3 and all(_even_profile(p)[1] == 0 and _even_profile(p)[0] > 10.0 for p in series.values())


def test_planarize_conjugated_nh2():
    np = pytest.importorskip("numpy")
    pytest.importorskip("pBabel")
    sys.path.insert(0, os.path.join(_REPO, "src", "graphics_engine", "src"))
    import pBabel  # noqa: F401
    from pMolecule import System
    from pScientific.Geometry3 import Coordinates3
    from pdynamo import dyff_prep as DP
    # . acetamide: C0(=O1)(N2)C3, N2 bonded to H4 H5 placed out of the plane
    system = System.FromAtoms([6, 8, 7, 6, 1, 1], bonds=[(0, 1), (0, 2), (0, 3), (2, 4), (2, 5)])
    xyz = [[0.0, 0.0, 0.0], [0.0, 1.23, 0.0], [1.15, -0.7, 0.0], [-1.3, -0.75, 0.0],
           [1.4, -1.2, 0.9], [2.0, -0.3, -0.6]]
    c = Coordinates3.WithExtent(6)
    for i, p in enumerate(xyz):
        for k in range(3):
            c[i, k] = p[k]
    system.coordinates3 = c
    assert DP.planarize_conjugated_nh2(system, built={4, 5}) == 1
    for h in (4, 5):
        assert abs(system.coordinates3[h, 2]) < 1e-6                  # . in the amide plane (z = 0)
    # . hydrogens from the input are not moved
    system.coordinates3[4, 2] = 0.5
    assert DP.planarize_conjugated_nh2(system, built=set()) == 0
    assert system.coordinates3[4, 2] == 0.5


def test_completed_side_chain_avoids_clash():
    # . 4oas: the Gln18 side chain is missing in the crystal; completed in the extended
    #   conformation it landed 1.5 A from the Arg97 guanidinium (user report: bad Arg)
    pytest.importorskip("pBabel")
    sys.path.insert(0, os.path.join(_REPO, "src", "graphics_engine", "src"))
    import contextlib, io
    import pBabel  # noqa: F401
    from pBabel import ImportSystem
    from pdynamo import dyff_prep as DP
    pdb = os.path.join(_REPO, "examples", "dyff_eh_system_preparation", "inputs", "4oas_raw.pdb")
    if not os.path.exists(pdb):
        pytest.skip("example input missing")
    saved = os.dup(1)
    with open(os.devnull, "w") as devnull:
        os.dup2(devnull.fileno(), 1)
        try:
            raw = ImportSystem(pdb, log=None)
            info = DP.analyze(raw); DP.remove_residues(info, list(info["covalent"])); DP.apply_site_protonation(info)
            system, report = DP.build(info, complete_truncated=True, ph=7.4)
        finally:
            os.dup2(saved, 1); os.close(saved)
    rot = {r: (before, after) for r, before, after in report["completed_rotamers"]}
    assert rot["A:GLN.18"][0] < 2.0 and rot["A:GLN.18"][1] > 3.5
    assert all(after >= before for before, after in rot.values())
