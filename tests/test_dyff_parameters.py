"""
Testes do backend da janela DYFF Parameters do Builder
(src/gui/windows/builder/dyff_parameters.py):

  1) os termos que o DYFF gera sao lidos do modelo construido, com a
     tensao (strain) de cada termo na geometria atual;
  2) um tipo atomico ERRADO (override manual absurdo) aparece como
     tensao alta (vermelho) -- e' o sinal de qualidade do DYFF;
  3) cargas parciais: o DYFF nao tem nenhuma (todas zero) -> atomos
     amarelos "no partial charges"; cargas definidas no Builder
     completam o conjunto.

Precisam do pDynamo3.
"""
import os
import sys

import numpy as np
import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_PDYNAMO = os.environ.get("PDYNAMO3_HOME", os.path.expanduser("~/programs/pDynamo3"))
if _PDYNAMO not in sys.path:
    sys.path.insert(0, _PDYNAMO)
pytest.importorskip("pBabel")
os.environ.setdefault("PDYNAMO3_PARAMETERS", os.path.join(_PDYNAMO, "parameters"))


@pytest.fixture
def adenosine(make_vismol_object):
    from pdynamo.opls import ligand
    atoms, bonds = ligand._parse_mol2(os.path.join(_TESTS_DIR, "data", "adenosine.mol2"))
    elements = [a[1] for a in atoms]
    kekule, _ = ligand._kekulize(elements, bonds)
    vobj = make_vismol_object(elements, [x for (i, j, o) in kekule for x in (i, j)])
    for k, a in enumerate(atoms):
        vobj.frames[0, k] = a[2]
        vobj.atoms[k].name = "%s%d" % (elements[k], k + 1)
    vobj.manual_bonds = {(min(i, j), max(i, j)) for (i, j, o) in kekule}
    vobj.manual_bond_orders = {(min(i, j), max(i, j)): o for (i, j, o) in kekule}
    vobj.manual_aromatic_bonds = set()
    return vobj


def test_terms_are_read_with_their_strain(adenosine):
    from gui.windows.builder import dyff_parameters as D
    tables = D.build_tables(adenosine)
    assert tables["error"] is None
    for term in ("bond", "angle", "dihedral", "outofplane"):
        assert tables[term], term
        for row in tables[term]:
            assert row["count"] == len(row["instances"]) and row["strain"] >= 0.0
    # . the crystal-like geometry agrees with DYFF: bonds are almost unstrained
    assert max(r["strain"] for r in tables["bond"]) < D.STRAIN_YELLOW
    # . DYFF has no partial charges: unstrained atoms are flagged for that
    assert not tables["has_charges"]
    assert any(r["quality"] == "no partial charges" and r["level"] == "yellow" for r in tables["atoms"])
    # . ring-closure angles are judged against the ring shape, not the natural angle
    assert any("ring" in r["quality"] for r in tables["angle"])


def test_wrong_atom_type_shows_up_as_strain(adenosine):
    from gui.windows.builder import dyff_parameters as D
    from gui.windows.builder.atom_types import set_manual_atom_type_override
    before = D.build_tables(adenosine)
    worst_before = max(r["strain"] for r in before["angle"])
    # . an sp3 ring carbon of the ribose typed as linear sp carbon
    carbon = next(r for r in before["atoms"] if r["type"] == "C:Tet")
    set_manual_atom_type_override(adenosine, carbon["atom_id"], "C:Lin")
    after = D.build_tables(adenosine)
    assert max(r["strain"] for r in after["angle"]) > max(worst_before, D.STRAIN_YELLOW)
    row = next(r for r in after["atoms"] if r["atom_id"] == carbon["atom_id"])
    assert row["level"] == "blue" and row["type"] == "C:Lin"


def test_partial_charges_complete_the_set(adenosine):
    from gui.windows.builder import dyff_parameters as D
    tables = D.build_tables(adenosine)
    assert not tables["has_charges"]
    D.set_partial_charge(adenosine, 0, -0.25)
    tables = D.build_tables(adenosine)
    assert tables["has_charges"] and tables["total_partial"] == pytest.approx(-0.25)


def test_dyff_parameter_edits(adenosine, tmp_path):
    from gui.windows.builder import dyff_parameters as D
    from gui.windows.builder.atom_ops import _snapshot_builder_state
    tables = D.build_tables(adenosine)
    row = max(tables["dihedral"], key=lambda r: r["count"])
    D.set_parameter_edit(adenosine, "dihedral", row["key"], [(5.0, 3)])
    edited = next(r for r in D.build_tables(adenosine)["dihedral"] if r["key"] == row["key"])
    assert edited["value"] == [(5.0, 3)] and edited["level"] == "blue" and edited["count"] == row["count"]
    assert _snapshot_builder_state(adenosine)["dyff_parameter_edits"]["dihedral"][row["key"]] == [[5.0, 3]]   # kept by Undo
    D.set_parameter_edit(adenosine, "dihedral", row["key"], None)                                          # back to the rule
    back = next(r for r in D.build_tables(adenosine)["dihedral"] if r["key"] == row["key"])
    assert back["value"] == row["value"] and not back["edited"]


def test_edited_dyff_model_energy_and_pickle(tmp_path):
    from pBabel import ImportSystem
    from pCore import Pickle, Unpickle
    from pMolecule.NBModel import NBModelCutOff
    from pdynamo.dyff_edits import dyff_model, MMModelDYFFEdited
    path = os.path.join(_TESTS_DIR, "data", "adenosine.mol2")
    plain = ImportSystem(path, log=None)
    plain.DefineMMModel(dyff_model("dyff-1.0"), log=None); plain.DefineNBModel(NBModelCutOff.WithDefaults())
    e_plain = plain.Energy(log=None)
    key = plain.mmState.mmTerms[[type(t).__name__ for t in plain.mmState.mmTerms].index("CosineDihedralContainer")].__getstate__()["parameterKeys"][0]
    edited = ImportSystem(path, log=None)
    edited.DefineMMModel(dyff_model("dyff-1.0", {"dihedral": {key: [[50.0, 1]]}}), log=None); edited.DefineNBModel(NBModelCutOff.WithDefaults())
    assert isinstance(edited.mmModel, MMModelDYFFEdited)
    e_edited = edited.Energy(log=None)
    assert abs(e_edited - e_plain) > 1.0
    Pickle(str(tmp_path / "s.pkl"), edited)
    again = Unpickle(str(tmp_path / "s.pkl"))
    assert again.Energy(log=None) == pytest.approx(e_edited)
