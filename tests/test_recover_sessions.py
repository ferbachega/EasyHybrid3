"""
Testes da janela "Unsaved Sessions" (src/gui/windows/setup/recover_sessions_window.py):
o preview de um .tmp/*.easy orfao precisa bater com o conteudo real e, principalmente,
ser rapido mesmo em arquivos grandes -- ver o modulo para o porque (pMolecule.System.
__setstate__ sempre refaz a percepcao de ordem de ligacao/aromaticidade do zero).
"""
import os
import pickle

import pytest

PDYNAMO_ENV = "/home/fernando/programs/pDynamo3/installation/shellScripts/environment_bash.com"


def _bala_system():
    import pBabel
    from pMolecule.MMModel import MMModelOPLS
    from pMolecule.NBModel import NBModelFull
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    system = pBabel.ImportSystem(os.path.join(repo, "examples", "tutorials", "01_bala", "bALA_not_opt.mol"))
    system.DefineMMModel(MMModelOPLS.WithParameterSet("protein"))
    system.DefineNBModel(NBModelFull.WithDefaults())
    return system


def test_describe_tmp_easy_file(tmp_path):
    pytest.importorskip("pBabel")
    if not os.path.exists(PDYNAMO_ENV):
        pytest.skip("pDynamo environment script not found")
    from gui.windows.setup.recover_sessions_window import describe_tmp_easy_file

    system = _bala_system()
    system.label, system.e_tag = "bALA test", "bALA"
    system.e_job_history = {0: {"simulation_type": "Geometry_Optimization", "status": "Finished"}}

    path = tmp_path / "orphan.easy"
    with open(path, "wb") as f:
        pickle.dump({"systems": [{"system": system, "vobjects": [{"name": "0_coords"}]}]}, f)

    text = describe_tmp_easy_file(str(path))
    assert "bALA test (bALA)" in text
    assert "22 atoms" in text and "force field: OPLS" in text
    assert "0_coords" in text
    assert "[0] Geometry_Optimization - Finished" in text


def test_describe_tmp_easy_file_unreadable(tmp_path):
    from gui.windows.setup.recover_sessions_window import describe_tmp_easy_file
    path = tmp_path / "not_a_pickle.easy"
    path.write_bytes(b"not a pickle at all")
    assert describe_tmp_easy_file(str(path)).startswith("Could not read this file:")


def test_fast_preview_unpickler_skips_connectivity(tmp_path):
    """ Regression guard for the speed fix: a System unpickled through
        _fast_preview_unpickler() must have no working connectivity
        (_connectivity is None, Connectivity.FromAtoms() never ran) --
        that is the entire point, and what makes previewing a big .tmp
        file fast instead of taking seconds. """
    pytest.importorskip("pBabel")
    if not os.path.exists(PDYNAMO_ENV):
        pytest.skip("pDynamo environment script not found")
    from gui.windows.setup.recover_sessions_window import _fast_preview_unpickler

    system = _bala_system()
    system.label, system.e_tag = "bALA test", "bALA"
    assert system.connectivity is not None                 # . the real System did perceive it

    path = tmp_path / "orphan.easy"
    with open(path, "wb") as f:
        pickle.dump({"systems": [{"system": system, "vobjects": []}]}, f)

    with open(path, "rb") as f:
        fast = _fast_preview_unpickler(f).load()
    fast_system = fast["systems"][0]["system"]
    assert fast_system.label == "bALA test" and len(fast_system.atoms) == 22
    assert fast_system._connectivity is None
