"""
Testes da fase 0 da preparacao OPLS (src/pdynamo/opls/):

  1) hydrogen_builder: prolina dentro da cadeia NAO recebe H amida
     (bug: o template da biblioteca e' o aminoacido livre, com H no N).
  2) prep: crambina so com atomos pesados -> 642 atomos, 3 pontes S-S,
     carga 0, 100 % tipada pelo conjunto OPLS "protein".
  3) prep: OXT ausente no C-terminal e' posicionado pela geometria.
  4) parameters: conversoes Tinker -> pDynamo (oplsaa08 reproduz quase
     todo o conjunto "protein" do pDynamo).
  5) parameters: o conjunto materializado para o sistema da EXATAMENTE a
     mesma energia que o conjunto "protein" original.

Precisam do pDynamo3 (PDYNAMO3_HOME ou ~/programs/pDynamo3); sem ele, os
testes sao pulados.
"""
import os
import sys
import shutil

import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_TESTS_DIR)
_PDYNAMO = os.environ.get("PDYNAMO3_HOME", os.path.expanduser("~/programs/pDynamo3"))
for _p in (_PDYNAMO, os.path.join(_REPO_ROOT, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

pytest.importorskip("pBabel")
if not os.getenv("PDYNAMO3_PARAMETERS"):
    os.environ["PDYNAMO3_PARAMETERS"] = os.path.join(_PDYNAMO, "parameters")

from pBabel import ImportSystem                       # noqa: E402
from pdynamo.opls import prep, parameters             # noqa: E402

CRAMBIN = os.path.join(_TESTS_DIR, "data", "crambin_heavy_atoms.pdb")


@pytest.fixture
def crambin():
    # . fresh system per test: hydrogen_builder.rebuild_system_with_added_hydrogens()
    #   reuses the Atom objects of the system it is given
    return ImportSystem(CRAMBIN, log=None)


def test_proline_gets_no_amide_hydrogen(crambin):
    from util import hydrogen_builder as hb
    results = hb.analyze_system(crambin)
    choices = {r["key"]: r["default_choice"] for r in results if r.get("is_standard_residue")}
    system = hb.rebuild_system_with_added_hydrogens(crambin, choices, random_seed=1)["system"]
    adjacent = system.connectivity.adjacentNodes
    for atom in system.atoms:
        if atom.label == "N" and atom.parent.label.startswith("PRO"):
            assert not any(n.atomicNumber == 1 for n in adjacent[atom]), atom.path
    assert len(system.atoms) == 642


def test_crambin_prepared_and_fully_typed(crambin):
    analysis = prep.analyze_system(crambin)
    assert len(analysis.disulfides) == 3
    assert not analysis.unsupported and not analysis.chain_breaks
    system, report = prep.build_prepared_system(analysis)
    assert report["atoms"] == 642
    assert report["formal_charge"] == 0
    assert not report["undefined_atoms"]
    typing = prep.typing_report(system)
    assert typing["ok"], typing["untyped_by_residue"]
    assert abs(typing["charge"]) < 1e-6


def test_missing_oxt_is_placed(tmp_path):
    path = tmp_path / "no_oxt.pdb"
    with open(CRAMBIN) as source, open(path, "w") as target:
        target.writelines(line for line in source if " OXT " not in line)
    analysis = prep.analyze_system(ImportSystem(str(path), log=None))
    assert analysis.summary()["added_atoms"] == {"A:ASN.46": ["OXT"]}
    system, report = prep.build_prepared_system(analysis)
    assert report["atoms"] == 642 and not report["undefined_atoms"]


def test_tinker_conversions_match_pdynamo():
    report = parameters.compare_with_pdynamo()
    bonds_eq, bonds_diff, _ = report["oplsaa08"]["bond"]
    angles_eq, angles_diff, _ = report["oplsaa08"]["angle"]
    assert bonds_eq >= 300 and bonds_diff <= 5
    assert angles_eq >= 850 and angles_diff <= 20
    oop_eq, oop_diff, _ = report["oplsaal"]["outofplane"]
    assert oop_eq >= 10 and oop_diff == 0


def test_materialized_set_reproduces_protein_energy(crambin, tmp_path):
    from pCore import Clone
    from pMolecule.MMModel import MMModelOPLS
    from pMolecule.NBModel import NBModelCutOff
    system, _ = prep.build_prepared_system(prep.analyze_system(crambin))
    types, _, untyped = prep.type_atoms(system)
    assert not untyped
    folder = str(tmp_path / "set")
    report = parameters.materialize_parameter_set(system, types, folder)
    assert not any(report["missing"].values())

    def energy(parameter_set):
        copy = Clone(system)
        copy.DefineMMModel(MMModelOPLS.WithParameterSet(parameter_set), log=None)
        copy.DefineNBModel(NBModelCutOff.WithDefaults())
        return copy.Energy(log=None)

    assert energy(os.path.abspath(folder)) == pytest.approx(energy("protein"), abs=1e-6)


# ---------------------------------------------------------------------------
#  Phase 2 (solvation) and phase 4 (ligands)
# ---------------------------------------------------------------------------
def test_solvated_cube_is_neutral_and_periodic(crambin, tmp_path):
    system, report = prep.prepare_opls_system(
        crambin, work_dir=str(tmp_path / "w"),
        solvation=dict(shape="cube", padding=6.0, concentration=0.15))
    assert report["ok"], report["typing"]["untyped_by_residue"]
    solvation = report["solvation"]
    assert solvation["waters"] > 500
    assert report["formal_charge"] == 0
    assert system.symmetry is not None
    assert abs(report["typing"]["charge"]) < 1e-6


def _have(tool):
    import shutil
    from pdynamo.opls import ligand
    return (ligand.find_xtb() if tool == "xtb" else shutil.which(tool)) is not None


@pytest.mark.skipif(not (_have("xtb") and _have("obabel")), reason="xtb and obabel needed")
def test_ligand_from_mol2_types_and_charges(tmp_path):
    system = ImportSystem(os.path.join(_TESTS_DIR, "data", "adenosine.mol2"), log=None)
    new_system, report = prep.prepare_opls_system(system, work_dir=str(tmp_path / "w"))
    assert report["ok"], report.get("typing")
    ligand = list(report["ligands"].values())[0]
    assert ligand["total_charge"] == 0 and abs(ligand["charge_sum"]) < 1e-6
    for opls_class in ("N*", "CK", "CB", "CQ", "NC", "N2"):          # purine nucleoside classes
        assert opls_class in ligand["classes"], opls_class
    assert not any(e.split()[0] in ("bond", "angle") for e in report["estimated_terms"])


@pytest.mark.skipif(not _have("obabel"), reason="obabel needed")
def test_cationic_aromatic_ring_gets_its_charge():
    from pdynamo.opls import ligand
    atoms = [(n, ligand.ELEMENT_Z[e], x) for (n, e, x, t) in ligand._parse_mol2(os.path.join(_TESTS_DIR, "data", "thiamine.mol2"))[0]]
    lig = ligand.ligand_chemistry("THI", atoms)
    assert lig.total_charge == +1                                      # thiazolium N+
    valence = {}
    for i, j, o in lig.bonds:
        valence[i] = valence.get(i, 0) + o; valence[j] = valence.get(j, 0) + o
    assert all(valence[k] == 4 for k, (n, z, x) in enumerate(lig.atoms) if z == 6)


# ---------------------------------------------------------------------------
#  Systems without residues (Builder, xyz), fragmented mol2, mol2 fallback,
#  quality report, class rules
# ---------------------------------------------------------------------------
def _builder_like_system():
    """ The calls of gui/windows/builder/empty_object._build_pdynamo_system_from_vismol_object:
        Connectivity + System.FromConnectivity, NO sequence (adenosine + one water). """
    import numpy as np
    from pMolecule import Atom, Bond, BondType, Connectivity, ConvertInputConnectivity, System
    from pScientific.Geometry3 import Coordinates3
    from pdynamo.opls import ligand
    atoms, bonds = ligand._parse_mol2(os.path.join(_TESTS_DIR, "data", "adenosine.mol2"))
    elements = [a[1] for a in atoms]
    kekule, _ = ligand._kekulize(elements, bonds)
    connectivity, xyz = Connectivity(), []
    for name, element, position, _ in atoms:
        connectivity.AddNode(Atom.WithOptions(atomicNumber=ligand.ELEMENT_Z[element], label=element))
        xyz.append(position)
    types = {1: BondType.Single, 2: BondType.Double, 3: BondType.Triple}
    for i, j, order in kekule:
        connectivity.AddEdge(Bond.WithNodes(connectivity.nodes[i], connectivity.nodes[j], type=types[order]))
    base = len(atoms)
    for z, position in ((8, (15.0, 0, 0)), (1, (15.96, 0, 0)), (1, (14.76, 0.93, 0))):
        connectivity.AddNode(Atom.WithOptions(atomicNumber=z, label="O" if z == 8 else "H"))
        xyz.append(np.array(position))
    for k in (1, 2):
        connectivity.AddEdge(Bond.WithNodes(connectivity.nodes[base], connectivity.nodes[base + k], type=BondType.Single))
    ConvertInputConnectivity(connectivity, {})
    system = System.FromConnectivity(connectivity=connectivity)
    coordinates = Coordinates3.WithExtent(len(xyz))
    for k, p in enumerate(xyz):
        coordinates[k, 0], coordinates[k, 1], coordinates[k, 2] = p
    system.coordinates3 = coordinates
    return system


def test_system_without_residues_is_split_into_molecules():
    system = _builder_like_system()
    analysis = prep.analyze_system(system)                    # used to raise "no residue information"
    kinds = [(r.name, r.kind) for r in analysis.residues]
    assert kinds == [("LIG", "unsupported"), ("HOH", "water")]


@pytest.mark.skipif(not (_have("xtb") and _have("obabel")), reason="xtb and obabel needed")
def test_builder_molecule_is_prepared(tmp_path):
    new_system, report = prep.prepare_opls_system(_builder_like_system(), work_dir=str(tmp_path / "w"))
    assert report["ok"] and report["atoms"] == 35
    assert "quality" in report


@pytest.mark.skipif(not (_have("xtb") and _have("obabel")), reason="xtb and obabel needed")
def test_mol2_fragments_are_merged_into_one_ligand(tmp_path):
    # . the pDynamo mol2 reader splits folate into "GLU"/"UNK" substructures
    system = ImportSystem(os.path.join(_TESTS_DIR, "data", "folic_acid.mol2"), log=None)
    analysis = prep.analyze_system(system)
    assert len(analysis.residues) == 1 and analysis.residues[0].kind == "unsupported"
    new_system, report = prep.prepare_opls_system(system, analysis=analysis, work_dir=str(tmp_path / "w"))
    assert report["ok"] and report["atoms"] == 51                # C19H19N7O6


@pytest.mark.skipif(not _have("obabel"), reason="obabel needed")
def test_mol2_fallback_reads_charged_aromatics():
    from pdynamo.opls import mol2_import
    system, used_fallback = mol2_import.import_system(os.path.join(_TESTS_DIR, "data", "thiamine.mol2"))
    assert used_fallback
    assert sum(a.formalCharge for a in system.atoms) == +1


@pytest.mark.skipif(not (_have("xtb") and _have("obabel")), reason="xtb and obabel needed")
def test_quality_report_flags_phosphate_torsions(tmp_path):
    system = ImportSystem(os.path.join(_TESTS_DIR, "data", "atp.mol2"), log=None)
    new_system, report = prep.prepare_opls_system(system, work_dir=str(tmp_path / "w"))
    assert report["ok"]
    quality = list(report["quality"].values())[0]
    assert quality["torsion_scans"]                             # P-O torsions have no OPLS parameters
    assert quality["levels"]["red"] > 0 and quality["charge_outliers"]
    assert os.path.exists(str(tmp_path / "w" / "quality.json"))


@pytest.mark.skipif(not (_have("xtb") and _have("obabel")), reason="xtb and obabel needed")
def test_imidazole_classes_match_opls():
    import subprocess, tempfile
    from pdynamo.opls import ligand
    scratch = tempfile.mkdtemp()
    subprocess.run(["obabel", "-:[nH]1cncc1", "--gen3d", "-h", "-O", os.path.join(scratch, "i.mol2")], capture_output=True)
    atoms, bonds = ligand._parse_mol2(os.path.join(scratch, "i.mol2"))
    elements = [a[1] for a in atoms]
    kekule, _ = ligand._kekulize(elements, bonds)
    lig = ligand.parametrize_residue("IMI", [("%s%d" % (e, k), ligand.ELEMENT_Z[e], a[2]) for k, (a, e) in enumerate(zip(atoms, elements))],
                                     bonds=kekule, cm5_scale=1.0)
    assert lig.classes[:5] == ["NA", "CR", "NB", "CV", "CW"]    # OPLS-AA imidazole N1 C2 N3 C4 C5


@pytest.mark.skipif(not (_have("xtb") and _have("obabel")), reason="xtb and obabel needed")
def test_ligand_charge_and_multiplicity_set_by_the_user(tmp_path):
    system = ImportSystem(os.path.join(_TESTS_DIR, "data", "adenosine.mol2"), log=None)
    analysis = prep.analyze_system(system)
    name = list(prep.ligand_preview(analysis))[0]
    analysis.ligand_options = {name: {"charge": 1, "multiplicity": 2}}        # radical cation
    new_system, report = prep.prepare_opls_system(system, analysis=analysis, work_dir=str(tmp_path / "w"),
                                                  solvation=dict(shape="cube", padding=6.0, concentration=0.0))
    assert report["ok"]
    ligand = report["ligands"][name]
    assert ligand["total_charge"] == 1 and ligand["multiplicity"] == 2
    assert report["solvation"]["ions"] == {"Cl-": 1}                           # neutralizes the user's +1
    assert report["total_charge"] == 0 and abs(report["typing"]["charge"]) < 1e-6


def test_incompatible_multiplicity_is_rejected():
    import numpy as np
    from pdynamo.opls import ligand
    lig = ligand.LigandParameters("MET")
    lig.atoms = [("C1", 6, np.zeros(3))] + [("H%d" % k, 1, np.zeros(3)) for k in range(1, 5)]
    lig.formal_charges = [0] * 5                                                # methane: 10 electrons
    with pytest.raises(ValueError):
        ligand.apply_user_charge(lig, total_charge=0, multiplicity=2)
    ligand.apply_user_charge(lig, total_charge=1, multiplicity=2)               # 9 electrons, doublet
    assert lig.total_charge == 1 and lig.n_electrons == 9


# ---------------------------------------------------------------------------
#  OPLS Parameters window backend (pdynamo/opls/editor.py)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not (_have("xtb") and _have("obabel")), reason="xtb and obabel needed")
def test_parameter_tables_and_editing(tmp_path):
    from pdynamo.opls import editor
    system = ImportSystem(os.path.join(_TESTS_DIR, "data", "atp.mol2"), log=None)
    new_system, report = prep.prepare_opls_system(system, work_dir=str(tmp_path / "w"))
    new_system.e_input_files = {"opls_parameters": report["parameter_folder"]}
    tables = editor.build_tables(new_system)
    red = [r for r in tables["dihedral"] if r["level"] == "red"]
    assert red and all(r["quality"] == "no torsion, rotatable bond" for r in red)
    q0, q1 = tables["charges"][0]["value"], tables["charges"][1]["value"]
    with pytest.raises(ValueError):                                     # non-integral total charge
        editor.apply_edits(new_system, {"charges": {0: q0 + 0.1}})
    editor.apply_edits(new_system, {"charges": {0: q0 + 0.1, 1: q1 - 0.1}})
    summary = editor.apply_edits(new_system, {"dihedral": {red[0]["key"]: [(0.5, 3, 0.0)]}})
    assert summary["rebuilt"] and new_system.nbModel is not None         # NB model put back
    after = editor.build_tables(new_system)
    assert [r["level"] for r in after["dihedral"] if r["key"] == red[0]["key"]] == ["blue"]
    assert abs(after["charges"][0]["value"] - (q0 + 0.1)) < 1e-6           # ligand charge survives the rebuild


def test_term_geometry_measures():
    from pdynamo.opls import editor
    assert editor.measure("bond", [(0, 0, 0), (1.5, 0, 0)]) == pytest.approx(1.5)
    assert editor.measure("angle", [(1, 0, 0), (0, 0, 0), (0, 1, 0)]) == pytest.approx(90.0)
    assert abs(editor.measure("dihedral", [(1, 0, 0), (0, 0, 0), (0, 1, 0), (0, 1, 1)])) == pytest.approx(90.0)
    assert editor.measure("outofplane", [(0, 0, 0.3), (1, 0, 0), (0, 1, 0), (-1, -1, 0)]) == pytest.approx(0.3)


def test_parameter_rows_carry_every_instance(crambin, tmp_path):
    from pdynamo.opls import editor
    system, report = prep.prepare_opls_system(crambin, work_dir=str(tmp_path / "w"))
    system.e_input_files = {"opls_parameters": report["parameter_folder"]}
    tables = editor.build_tables(system)
    for term, size in (("bond", 2), ("angle", 3), ("dihedral", 4), ("outofplane", 4)):
        for row in tables[term]:
            assert len(row["instances"]) == row["count"] and all(len(t) == size for t in row["instances"])
    assert sum(len(r["instances"]) for r in tables["lj"]) == len(system.atoms)
