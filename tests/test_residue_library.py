"""
Biblioteca de residuos do DYFF (src/pdynamo/residue_library.py,
src/pdynamo/data/dyff_residue_library.json):

  1) a biblioteca tem os 71 modelos completos (26 internos, 44 terminais N/C,
     CYX), cada residuo com carga inteira;
  2) o reconhecimento pelo grafo escolhe a variante de protonacao certa
     (HID x HIE x HIP; ASP x ASH) a partir dos atomos presentes;
  3) um Cys em ponte dissulfeto NAO e' confundido com o tiolato CYM: casa com
     CYX, com a juncao x:S no S do parceiro; o tiolato livre casa com CYM;
  4) crambina: todos os 46 residuos cobertos (Thr N-terminal, Asn C-terminal,
     6 CYX), edicoes so com indices validos.
"""
import os
import sys

import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_TESTS_DIR)
_PDYNAMO = os.environ.get("PDYNAMO3_HOME", os.path.expanduser("~/programs/pDynamo3"))
for p in (_PDYNAMO, os.path.join(_REPO, "src")):
    if p not in sys.path:
        sys.path.insert(0, p)
os.environ.setdefault("PDYNAMO3_PARAMETERS", os.path.join(_PDYNAMO, "parameters"))
pytest.importorskip("networkx")

from pdynamo import residue_library as RL     # noqa: E402


def _residue_from_template(t, offset=0):
    """ A residue (atoms offset..offset+n-1) with exactly the template's graph. """
    symbols = [a["element"] for a in t["atoms"]]
    neighbours = [set(offset + n for n in a["neighbors"]) for a in t["atoms"]]
    return symbols, neighbours


def test_library_is_complete_with_integer_residue_charges():
    lib = RL.load_library()
    assert len(lib["templates"]) == 71
    assert sum(1 for n in lib["templates"] if len(n) == 4 and n[0] in "NC") == 44
    for name, t in lib["templates"].items():
        total = sum(a["charge"] for a in t["atoms"])
        assert abs(total - t["formal_charge"]) < 1e-3, name
        assert t["bond_terms"] and t["angle_terms"] and t["dihedral_bonds"], name


def test_protonation_variant_is_recognised_from_the_graph():
    lib = RL.load_library()
    for variant in ("HID", "HIE", "HIP", "ASP", "ASH"):
        symbols, neighbours = _residue_from_template(lib["templates"][variant])
        family = lib["templates"][variant]["family"]
        matches, unmatched = RL.match_residues([("A:{}.1".format(family), family, list(range(len(symbols))))],
                                               symbols, neighbours, lib)
        assert not unmatched and matches[0]["template"] == variant, variant


def test_disulfide_cys_is_not_a_thiolate():
    lib = RL.load_library()
    symbols, neighbours = _residue_from_template(lib["templates"]["CYM"])
    sg = next(k for k, a in enumerate(lib["templates"]["CYM"]["atoms"]) if a["element"] == "S")
    # . the free thiolate matches CYM ...
    m, u = RL.match_residues([("A:CYS.1", "CYS", list(range(len(symbols))))], symbols, neighbours, lib)
    assert m and m[0]["template"] == "CYM"
    # . ... but the same S bonded to another residue's S does not
    symbols = symbols + ["S"]
    neighbours = [set(n) for n in neighbours] + [{sg}]
    neighbours[sg].add(len(symbols) - 1)
    m, u = RL.match_residues([("A:CYS.1", "CYS", list(range(len(symbols) - 1)))], symbols, neighbours, lib)
    assert m and m[0]["template"] == "CYX" and m[0]["junction"]["x:S"] == len(symbols) - 1


def test_crambin_plan():
    pytest.importorskip("pBabel")
    import contextlib, io, tempfile
    from pBabel import ImportSystem
    from pdynamo.opls import prep
    with contextlib.redirect_stdout(io.StringIO()):
        heavy = ImportSystem(os.path.join(_TESTS_DIR, "data", "crambin_heavy_atoms.pdb"), log=None)
        system, report = prep.prepare_opls_system(heavy, work_dir=tempfile.mkdtemp())
    p = RL.plan(system)
    s = p["summary"]
    assert s["residues"] == 46 and s["matched"] == 46 and not p["rest_atoms"]
    assert s["templates"].get("CYX") == 6 and s["templates"].get("NTHR") == 1 and s["templates"].get("CASN") == 1
    n = len(system.atoms)
    for kind in ("bond_terms", "angle_terms", "dihedral_bonds"):
        for key in p["edits"][kind]:
            assert all(0 <= int(i) < n for i in key.split("-"))
    assert abs(sum(p["charges"].values()) - round(sum(p["charges"].values()))) < 1e-3


def test_water_and_ions_get_fixed_charges():
    """ Water (O with exactly two H) and monatomic ions are recognised by the graph and
        get TIP3P / formal charges; a hydroxide or an alcohol O is not water. """
    symbols = ["O", "H", "H", "Na", "Cl", "O", "H", "C", "O", "H"]
    neighbours = [{1, 2}, {0}, {0}, set(), set(), {6}, {5}, {8}, {7, 9}, {8}]
    water, ions, charges = RL.solvent_and_ions(symbols, neighbours, range(len(symbols)))
    assert water == [0, 1, 2] and ions == [3, 4]
    assert charges[0] == pytest.approx(-0.834) and charges[1] == pytest.approx(0.417)
    assert charges[3] == 1.0 and charges[4] == -1.0
    assert 5 not in charges and 8 not in charges
