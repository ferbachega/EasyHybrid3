"""
Ajuste de diedros OPLS contra uma varredura xTB
(src/pdynamo/opls/torsion_fit.py + src/pdynamo/torsion_fit.py):

  1) a energia MM com o parametro trocado (fragmento podado - termo atual +
     serie ajustada) e' igual, a menos de uma constante, a energia do modelo
     MM RECONSTRUIDO a partir dos arquivos editados (editor.apply_edits);
  2) conversao serie de cossenos <-> forma OPLS k ( 1 + cos ( n phi - delta ) );
  3) varredura real (xtb) de uma ligacao giravel da adenosina: o ajuste
     reduz o erro e o parametro aplicado aparece como editado (azul).

Precisam do pDynamo3, obabel e xtb.
"""
import math
import os
import sys

import numpy as np
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

from pdynamo import torsion_fit as F                  # noqa: E402
from pdynamo.opls import torsion_fit as OT            # noqa: E402


def _have_tools():
    import shutil
    from pdynamo.opls import ligand
    return ligand.find_xtb() is not None and shutil.which("obabel") is not None


needs_tools = pytest.mark.skipif(not _have_tools(), reason="xtb and obabel needed")


@pytest.fixture(scope="module")
def adenosine_opls(tmp_path_factory):
    from pBabel import ImportSystem
    from pdynamo.opls import prep, editor
    work = tmp_path_factory.mktemp("ado")
    system = ImportSystem(os.path.join(_TESTS_DIR, "data", "adenosine.mol2"), log=None)
    new_system, report = prep.prepare_opls_system(system, work_dir=str(work / "w"))
    new_system.e_input_files = {"opls_parameters": report["parameter_folder"]}
    return new_system


def _rotatable_row(system, tables):
    """ a dihedral row with a term around a rotatable (acyclic, non-terminal) bond. """
    n = len(system.atoms)
    bonds = list(system.connectivity.bondIndices)
    adjacency = F.adjacency_from_bonds(n, bonds)
    for row in tables["dihedral"]:
        for q in row["instances"]:
            j, k = q[1], q[2]
            if len(adjacency[j]) < 2 or len(adjacency[k]) < 2:
                continue
            try:
                F.moving_side(adjacency, j, k)
            except ValueError:
                continue
            if system.atoms[q[0]].atomicNumber > 1 and system.atoms[q[3]].atomicNumber > 1:
                return row, q
    raise AssertionError("no rotatable dihedral")


def test_fourier_conversion():
    pairs = [(4.184, 1), (-8.368, 2), (2.092, 3)]
    value = OT.pairs_to_fourier(pairs)
    assert value == [(1.0, 1, 0.0), (2.0, 2, 180.0), (0.5, 3, 0.0)]
    for phi in (-150.0, -30.0, 0.0, 75.0, 180.0):
        series = sum(c * math.cos(n * math.radians(phi)) for c, n in pairs)
        constant = OT.KCAL_TO_KJ * sum(k for k, n, p in value)
        assert abs(OT.fourier_energy(value, phi) - (series + constant)) < 1e-9


@needs_tools
def test_replaced_parameter_energy_matches_a_real_rebuild(adenosine_opls, tmp_path):
    import shutil
    from pdynamo.opls import editor
    system = adenosine_opls
    tables = editor.build_tables(system)
    row, quad = _rotatable_row(system, tables)
    target = OT.OPLSTorsionTarget(system, row, quad)
    j, k = target.quad[1], target.quad[2]
    side = F.moving_side(F.adjacency_from_bonds(len(target.symbols), target.bonds), j, k)
    geometries = [F.rotate_about_bond(target.coordinates, j, k, side, a) for a in range(0, 360, 45)]
    pairs = [(3.0, 1), (-5.0, 2), (1.2, 3)]
    predicted = np.asarray(target.energies_for(pairs, geometries))
    # . the same edit through the files: a copy of the system's parameter set
    import copy
    folder = tables["folder"]
    copy_folder = str(tmp_path / "params")
    shutil.copytree(folder, copy_folder)
    clone = copy.deepcopy(system)
    clone.e_input_files = {"opls_parameters": copy_folder}
    editor.apply_edits(clone, {"dihedral": {target.key: OT.pairs_to_fourier(pairs)}})
    rebuilt = OT.OPLSTorsionTarget(clone, row, quad)
    real = np.asarray(F.mm_energies(rebuilt.fragment, geometries))
    difference = real - predicted
    # . constant offset only (the k of every term); k is written with 6 decimals
    assert np.ptp(difference) < 1e-4, difference


@needs_tools
def test_adenosine_torsion_fit_and_apply(adenosine_opls):
    from pdynamo.opls import editor
    system = adenosine_opls
    tables = editor.build_tables(system)
    row, quad = _rotatable_row(system, tables)
    target = OT.OPLSTorsionTarget(system, row, quad)
    runner = F.XTBRunner(gfn=2, charge=target.default_charge)
    scan = F.run_scan(target.symbols, target.coordinates, target.bonds, target.quad, step=30.0, relaxed=True, runner=runner)
    result = F.fit_dihedral_key(scan, target.key, target.instances, lambda p: target.energies_for(p, scan["geometries"]),
                                current_text=target.current_text, format_value=target.format_value, max_period=3,
                                force_field="OPLS")
    print(F.report_text(result))
    assert result["rmse_after"] < result["rmse_before"]
    message = target.apply(result["pairs"])
    after = editor.build_tables(system)
    assert [r["level"] for r in after["dihedral"] if tuple(r["key"]) == target.key] == ["blue"]
    assert system.nbModel is not None
