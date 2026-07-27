"""
Fixtures compartilhadas pelos testes de EasyHybrid3.

O objetivo aqui e' permitir testar a LOGICA DE DADOS do VismolObject
(bonds, percepcao de ordem de ligacao, Dynamic Bonds, etc) sem precisar
de GTK, OpenGL, ou um display real -- o que torna esses testes rodaveis
em CI (GitHub Actions) num runner comum, sem Xvfb nem nada parecido.

Isso e' possivel porque vismol.core.vismol_object e' importavel sozinho
desde que:
  1) as extensoes Cython do submodulo (src/graphics_engine) ja tenham
     sido compiladas -- ver src/graphics_engine/install.sh
     (`python3 setup.py build_ext --inplace`), que e' o mesmo comando
     que install.py ja usa em uma instalacao normal; e
  2) numpy, Cython, freetype-py e PyOpenGL estejam instalados (ja estao
     em requirements.txt) -- PyGObject/GTK3 NAO sao necessarios so' para
     isso, pois vismol_object.py nao importa gi diretamente.

VismolObject exige um vm_session/vm_config com alguns atributos
minimos (vm_config.gl_parameters, vm_config.representations_available,
periodic_table, vm_glcore) -- _FakeVmConfig/_FakeVmSession abaixo
fornecem exatamente o minimo necessario, sem depender de uma sessao
grafica completa.
"""
import os
import sys

import numpy as np
import pytest

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_TESTS_DIR)
for _p in (
    os.path.join(_REPO_ROOT, "src"),
    os.path.join(_REPO_ROOT, "src", "graphics_engine", "src"),
):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from vismol.core.vismol_object import VismolObject  # noqa: E402
from vismol.model.atom import Atom  # noqa: E402
from vismol.utils.elements import PeriodicTable  # noqa: E402


class _FakeVmConfig:
    """ So' os atributos que VismolObject.__init__ realmente le. """
    gl_parameters = {"multiple_bonds": True, "bond_tolerance": 1.4}
    representations_available = []


class _FakeVmSession:
    """ Substituto minimo de VismolSession -- fornece so' o que
    VismolObject/atom_ops.py precisam (vm_config, periodic_table,
    vm_glcore, get_frame/set_frame para os testes de Dynamic Bonds),
    sem nenhuma dependencia de GTK/OpenGL/display. """

    def __init__(self):
        self.vm_config = _FakeVmConfig()
        self.periodic_table = PeriodicTable()
        self.vm_glcore = None  # ausencia proposital -- sem contexto GL real
        self.atom_id_counter = 0
        self._frame = 0

    def get_frame(self):
        return self._frame

    def set_frame(self, f):
        self._frame = int(f)


@pytest.fixture
def vm_session():
    return _FakeVmSession()


@pytest.fixture
def make_vismol_object(vm_session):
    """ Factory fixture: make_vismol_object(symbols, bonds, n_frames=1)
    devolve um VismolObject pronto (atoms/frames/index_bonds/bonds ja
    montados, equivalente a uma estrutura recem-carregada de um arquivo
    real -- NAO construida atomo-a-atomo via add_atom()/add_bond(), que
    e' o cenario onde encontramos os bugs de "add_bond apaga as outras
    ligacoes" -- ver test_atom_ops_bonds.py).

    symbols : lista de simbolos quimicos, 1 por atomo (ex.: ['C','C','H']).
    bonds   : lista/array achatada de pares [i0,j0, i1,j1, ...].
    n_frames: numero de frames da "trajetoria" (frames sao so' zeros --
              estes testes nao dependem de coordenadas reais).

    create_representation() e' substituida por um stub no-op (nao
    depende de OpenGL) -- os testes verificam o estado de
    vismol_object.bonds/index_bonds/dynamic_bonds diretamente, que e'
    onde os bugs que motivaram esta suite realmente vivem. """
    def _make(symbols, bonds, n_frames=1):
        vobj = VismolObject(vismol_session=vm_session, name="test", active=True)
        vobj.vm_session = vm_session

        vobj.atoms = {}
        for i, sym in enumerate(symbols):
            vobj.atoms[i] = Atom(vismol_object=vobj, atom_id=i, name=sym,
                                  symbol=sym, index=i)
        vobj.frames = np.zeros((n_frames, len(symbols), 3), dtype=np.float32)

        vobj.index_bonds = list(bonds)
        vobj.bonds = None
        vobj._bonds_from_pair_of_indexes_list()

        vobj.create_representation = lambda rep_type=None, indexes=None: None
        vobj.core_representations = {"picking_dots": None, "picking_text": None}
        return vobj
    return _make


def pairs_in(flat_bonds):
    """ Helper: converte um array/lista achatada [i0,j0,i1,j1,...] num
    set de pares normalizados {(min,max), ...}, para comparar
    conectividade sem se importar com ordem/orientacao. """
    flat = list(np.asarray(flat_bonds).ravel().tolist()) if not isinstance(flat_bonds, list) else flat_bonds
    return {tuple(sorted((flat[k], flat[k + 1]))) for k in range(0, len(flat), 2)}
