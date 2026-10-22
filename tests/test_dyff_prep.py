"""
Preparo direto para DYFF + biblioteca de residuos, sem OPLS (src/pdynamo/dyff_prep.py),
com a base de metaloproteinas (tests/data/metalloproteins).

  1) 1ZNF (Zn Cys2His2): sitio detectado com os 4 doadores; as capas ACE/NH2,
     ligadas covalentemente, sao detectadas e podem ser removidas; o sistema tem as
     ligacoes Zn-doador, as Cys sao tiolato (sem HG), as His sem H no N ligado;
  2) residuo truncado (ubiquitina com a cadeia lateral de uma Lys cortada): os
     atomos faltantes sao recolocados com comprimentos de ligacao razoaveis;
  3) controle sem metal (1UBQ): nenhum sitio.
"""
import contextlib
import io
import os
import sys

import numpy as np
import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_TESTS_DIR)
for p in (os.path.join(_REPO, "src"), os.path.join(_REPO, "src", "graphics_engine", "src")):
    if p not in sys.path:
        sys.path.insert(0, p)
DATA = os.path.join(_TESTS_DIR, "data", "metalloproteins")
pytest.importorskip("pBabel")


def _import(path):
    import pBabel  # noqa: F401
    from pBabel import ImportSystem
    with contextlib.redirect_stdout(io.StringIO()):
        return ImportSystem(path, log=None)


def _neighbours(system):
    nb = {}
    for i, j in system.connectivity.bondIndices:
        nb.setdefault(i, []).append(j)
        nb.setdefault(j, []).append(i)
    return nb


def test_zinc_finger_site():
    from pdynamo import dyff_prep as DP
    info = DP.analyze(_import(os.path.join(DATA, "1ZNF.pdb")))
    residues = info["analysis"].residues
    assert len(info["metals"]) == 1
    site = info["metals"][0]
    assert site["element"] == "Zn" and site["model"] == "bonded"
    donors = sorted((residues[j].library_name, n) for j, n, d in site["donors"])
    assert donors == [("CYS", "SG"), ("CYS", "SG"), ("HIS", "NE2"), ("HIS", "NE2")]
    assert sorted(residues[i].name for i in info["covalent"]) == ["ACE", "NH2"]
    with pytest.raises(ValueError):
        DP.build(info)                                         # caps not removed
    DP.remove_residues(info, info["covalent"])
    DP.apply_site_protonation(info)
    with contextlib.redirect_stdout(io.StringIO()):
        system, report = DP.build(info)
    assert not report["undefined"]
    s = report["metal_sites"][0]
    nb = _neighbours(system)
    zn = s["metal_index"]
    assert sorted(system.atoms[k].label for k in nb[zn]) == ["NE2", "NE2", "SG", "SG"]
    for k in nb[zn]:
        hydrogens = [h for h in nb[k] if system.atoms[h].atomicNumber == 1]
        assert not hydrogens                                   # thiolate S, unprotonated N
    assert DP.site_groups(system, report)[zn]["multiplicity"] == 1


def test_truncated_residue_is_completed(tmp_path):
    from pdynamo import dyff_prep as DP
    lines = open(os.path.join(DATA, "1UBQ.pdb")).read().splitlines()
    # . Lys 6: keep N CA C O CB, drop CG CD CE NZ; drop the waters
    cut = [l for l in lines if l.startswith(("ATOM", "HETATM")) and not (l[17:20] == "LYS" and l[22:26].strip() == "6" and l[12:16].strip() in ("CG", "CD", "CE", "NZ"))
           and l[17:20] != "HOH"]
    path = tmp_path / "ubq_cut.pdb"
    path.write_text("\n".join(cut) + "\nEND\n")
    info = DP.analyze(_import(str(path)))
    assert [m for _, l, m in info["truncated"] if l.endswith("LYS.6")] == [["CG", "CD", "CE", "NZ"]]
    assert not info["metals"]
    with contextlib.redirect_stdout(io.StringIO()):
        system, report = DP.build(info)
    assert len(report["completed"]) == 4 and not report["undefined"]
    x = np.array([[system.coordinates3[a, k] for k in range(3)] for a in range(len(system.atoms))])
    idx = {a.label: a.index for a in system.atoms if a.path.split(":")[1] == "LYS.6"}
    for a, b in (("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")):
        assert 1.3 < np.linalg.norm(x[idx[a]] - x[idx[b]]) < 1.7
    # . nothing placed on top of another atom
    others = [i for i in range(len(x)) if i not in idx.values()]
    for name in ("CG", "CD", "CE", "NZ"):
        assert np.min(np.linalg.norm(x[others] - x[idx[name]], axis=1)) > 0.9


def test_untyped_metal_gets_a_type():
    """ Free Ca2+ (calmodulin, ionic sites) has no DYFF pattern: the preparation's type
        overrides make MMModelDYFFEdited type it ("Ca:Oct:+2"); without them typing fails. """
    from pdynamo import dyff_prep as DP
    from pdynamo.dyff_edits import dyff_model
    from pMolecule.MMModel.MMModelError import MMModelError
    info = DP.analyze(_import(os.path.join(DATA, "1CLL.pdb")))
    assert {s["element"] for s in info["metals"]} == {"Ca"} and all(s["model"] == "ionic" for s in info["metals"])
    DP.remove_residues(info, info["covalent"])
    with contextlib.redirect_stdout(io.StringIO()):
        system, report = DP.build(info)
    assert set(report["type_overrides"].values()) == {"Ca:Oct:+2"} and len(report["type_overrides"]) == 4
    assert report["repaired_hydrogens"] > 0                     # distorted crystal geometry repaired
    import copy
    plain = copy.deepcopy(system)
    with pytest.raises(MMModelError), contextlib.redirect_stdout(io.StringIO()):
        plain.DefineMMModel(dyff_model("dyff-1.0", None), log=None)
    with contextlib.redirect_stdout(io.StringIO()):
        system.DefineMMModel(dyff_model("dyff-1.0", None, report["type_overrides"]), log=None)


def test_altloc_pdb_through_the_import_path(tmp_path):
    """ 3KS3 (carbonic anhydrase II: HIS 4, waters 353/649 with ALTLOC A/B) imported the
        way File > Open does it (dedupe_pdb_atom_names first): the build gives every atom
        coordinates. Before the altLoc fix all 479 waters came out at the undefined
        sentinel (1e4 A), and the view centred on empty space. """
    import shutil
    from pdynamo import dyff_prep as DP
    from util.pdb_tools import dedupe_pdb_atom_names
    path = tmp_path / "3KS3.pdb"
    shutil.copy(os.path.join(DATA, "3KS3.pdb"), path)
    dedupe_pdb_atom_names(str(path))
    raw = _import(str(path))
    assert [a.label for a in raw.atoms if a.parent.label == "HOH.353"] == ["O"]
    info = DP.analyze(raw)
    DP.remove_residues(info, info["covalent"])
    DP.apply_site_protonation(info)
    with contextlib.redirect_stdout(io.StringIO()):
        system, report = DP.build(info)
    assert not report["undefined"]
    x = np.array([[system.coordinates3[a, k] for k in range(3)] for a in range(len(system.atoms))])
    assert np.abs(x).max() < 100.0


@pytest.mark.parametrize("pdb", ["1ZNF", "3KS3", "1CLL"])
def test_sites_from_structure_match_the_preparation(pdb):
    """ dyff_prep.sites_from_system (used when DYFF + residue library is applied again to a
        system edited in the Builder) reads the metal types and site states from the
        structure alone: on a freshly prepared system they equal the preparation's own. """
    from pdynamo import dyff_prep as DP
    info = DP.analyze(_import(os.path.join(DATA, pdb + ".pdb")))
    DP.remove_residues(info, info["covalent"])
    DP.apply_site_protonation(info)
    with contextlib.redirect_stdout(io.StringIO()):
        system, report = DP.build(info)
    overrides, states = DP.sites_from_system(system)
    assert {k: overrides[k] for k in report["type_overrides"]} == report["type_overrides"]
    assert states == DP.site_groups(system, report)
    # . the donor fallbacks are used only for untyped atoms: same types as the preparation
    import copy
    from pdynamo.dyff_edits import dyff_model
    types = []
    for o in (report["type_overrides"], overrides):
        s = copy.deepcopy(system)
        with contextlib.redirect_stdout(io.StringIO()):
            s.DefineMMModel(dyff_model("dyff-1.0", None, o), log=None)
        types.append(list(s.mmState.atomTypeIndices))
    assert types[0] == types[1]
    path = system.atoms[next(iter(states))].path
    assert DP.sites_from_system(system, {path: 3})[1][next(iter(states))]["multiplicity"] == 3


def test_carboxylate_donor_is_typed():
    """ Thermolysin (8TLN): Zn bonded to Glu166 OE2. The carboxylate O- bonded to the metal has
        no DYFF pattern; dyff_prep.type_overrides gives it a fallback type (used only for
        untyped atoms). Before the fix DYFF typing failed for every such site. """
    from pdynamo import dyff_prep as DP
    from pdynamo.dyff_edits import dyff_model
    info = DP.analyze(_import(os.path.join(DATA, "8TLN.pdb")))
    DP.remove_residues(info, info["covalent"])
    DP.apply_site_protonation(info)
    with contextlib.redirect_stdout(io.StringIO()):
        system, report = DP.build(info)
    zn = [s for s in report["metal_sites"] if s["element"] == "Zn" and s["model"] == "bonded"]
    assert zn
    donors = [system.atoms[d["index"]] for d in zn[0]["donors"] if d["index"] is not None]
    assert any(a.label.startswith("OE") for a in donors)
    with contextlib.redirect_stdout(io.StringIO()):
        system.DefineMMModel(dyff_model("dyff-1.0", None, report["type_overrides"]), log=None)
