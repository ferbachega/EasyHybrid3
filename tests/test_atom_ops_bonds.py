"""
Testes para as funcoes de edicao de ligacao em
gui/windows/builder/atom_ops.py: set_bond_order, unset_bond,
set_dynamic_bond_order, unset_dynamic_bond, resolve_frame_arg.

Cobrem bugs reais encontrados via teste manual/uso real do terminal
(ver o historico de commits em atom_ops.py e vismol_object.py para o
contexto completo):

  1) add_bond()/remove_bond() (as funcoes originais do Builder) reescrevem
     toda a conectividade a partir de vismol_object.manual_bonds -- que
     so' contem os pares explicitamente adicionados via add_bond(). Numa
     estrutura CARREGADA de arquivo (nao construida atomo-a-atomo no
     Builder), isso apaga todas as OUTRAS ligacoes da estrutura, que
     nunca estiveram em manual_bonds. set_bond_order()/unset_bond()
     foram criadas para corrigir isso, editando vismol_object.bonds
     diretamente.
  2) A representacao "sticks" (a unica que desenha duplas/triplas de
     forma diferenciada) nao era recriada por essas funcoes -- so'
     "lines"/"nonbonded" eram. _refresh_bond_dependent_representations()
     corrige isso.
  3) Dynamic Bonds (dynamic_bonds[f], por frame): set_dynamic_bond_order/
     unset_dynamic_bond precisam editar SO' os frames pedidos, sem afetar
     os demais, e aplicar overrides manuais de ordem por cima da
     percepcao automatica.
"""
import os
import sys

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_TESTS_DIR)
for _p in (
    os.path.join(_REPO_ROOT, "src"),
    os.path.join(_REPO_ROOT, "src", "graphics_engine", "src"),
):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from gui.windows.builder.atom_ops import (  # noqa: E402
    set_bond_order,
    unset_bond,
    resolve_frame_arg,
    set_dynamic_bond_order,
    unset_dynamic_bond,
)

from conftest import pairs_in  # noqa: E402


# ---------------------------------------------------------------------
# Helpers de fixture: uma estrutura "carregada de arquivo" -- benzeno
# (6C + 6H), SEM passar por add_atom()/add_bond() -- ou seja,
# vismol_object.manual_bonds nunca existiu, exatamente como um objeto
# vindo de um PDB de verdade.
# ---------------------------------------------------------------------

_BENZENE_SYMBOLS = ["C"] * 6 + ["H"] * 6
_BENZENE_BONDS = [
    0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0,
    0, 6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11,
]


def test_set_bond_order_preserves_every_other_bond(make_vismol_object):
    """ Bug principal: editar UMA ligacao numa estrutura carregada nao
    pode apagar as outras 11. """
    vobj = make_vismol_object(_BENZENE_SYMBOLS, _BENZENE_BONDS)
    n_bonds_before = len(vobj.bonds)
    assert n_bonds_before == 12

    set_bond_order(vobj, 0, 1, bond_order=2)

    assert len(vobj.bonds) == n_bonds_before, "nenhuma outra ligacao deveria ter sido perdida"
    assert vobj.bonds[(0, 1)].bond_order == 2


def test_unset_bond_removes_only_the_requested_pair(make_vismol_object):
    vobj = make_vismol_object(_BENZENE_SYMBOLS, _BENZENE_BONDS)
    n_bonds_before = len(vobj.bonds)

    removed = unset_bond(vobj, 2, 3)

    assert removed is True
    assert len(vobj.bonds) == n_bonds_before - 1
    assert (2, 3) not in vobj.bonds
    # as outras 10 ligacoes do anel/hidrogenios continuam intactas
    remaining_pairs = set(vobj.bonds.keys())
    expected = pairs_in(_BENZENE_BONDS) - {(2, 3)}
    assert remaining_pairs == expected


def test_unset_bond_on_nonexistent_pair_is_a_no_op(make_vismol_object):
    vobj = make_vismol_object(_BENZENE_SYMBOLS, _BENZENE_BONDS)
    n_bonds_before = len(vobj.bonds)

    removed = unset_bond(vobj, 0, 3)  # atomos opostos do anel, nunca ligados

    assert removed is False
    assert len(vobj.bonds) == n_bonds_before


def test_set_bond_order_creates_new_bond_without_disturbing_others(make_vismol_object):
    vobj = make_vismol_object(_BENZENE_SYMBOLS, _BENZENE_BONDS)
    n_bonds_before = len(vobj.bonds)

    created = set_bond_order(vobj, 0, 3, bond_order=1)  # atomos opostos, nunca ligados

    assert created is True
    assert len(vobj.bonds) == n_bonds_before + 1
    assert (0, 3) in vobj.bonds
    for pair in pairs_in(_BENZENE_BONDS):
        assert pair in vobj.bonds, f"ligacao original {pair} nao deveria ter sumido"


def test_set_bond_order_rejects_invalid_order(make_vismol_object):
    vobj = make_vismol_object(_BENZENE_SYMBOLS, _BENZENE_BONDS)
    import pytest
    with pytest.raises(ValueError):
        set_bond_order(vobj, 0, 1, bond_order=4)


def test_set_bond_order_rejects_self_bond(make_vismol_object):
    vobj = make_vismol_object(_BENZENE_SYMBOLS, _BENZENE_BONDS)
    import pytest
    with pytest.raises(ValueError):
        set_bond_order(vobj, 0, 0, bond_order=1)


# ---------------------------------------------------------------------
# Representacao "sticks" deve ser atualizada (bug: so' lines/nonbonded
# eram recriadas antes)
# ---------------------------------------------------------------------

class _FakeStickRepresentation:
    """ Dublê minimo de SticksRepresentation, so' com o que os testes
    precisam observar (indexes/active), sem nenhuma dependencia de
    OpenGL. """
    def __init__(self, indexes, active=True):
        self.indexes = list(indexes)
        self.active = active
    is_dynamic = False


def test_sticks_representation_is_refreshed_on_unbond(make_vismol_object):
    vobj = make_vismol_object(["C"] * 4, [0, 1, 1, 2, 2, 3])

    vobj.representations["sticks"] = _FakeStickRepresentation(indexes=[0, 1, 1, 2, 2, 3], active=True)
    vobj.representations["lines"] = _FakeStickRepresentation(indexes=[0, 1, 1, 2, 2, 3], active=False)

    def _fake_create_representation(rep_type=None, indexes=None):
        if rep_type == "sticks":
            flat = [v for pair in sorted(pairs_in(vobj.index_bonds)) for v in pair]
            vobj.representations["sticks"] = _FakeStickRepresentation(indexes=flat, active=True)
        elif rep_type == "lines":
            flat = [v for pair in sorted(pairs_in(vobj.index_bonds)) for v in pair]
            vobj.representations["lines"] = _FakeStickRepresentation(indexes=flat, active=True)
        # "nonbonded": no-op para este teste

    vobj.create_representation = _fake_create_representation

    unset_bond(vobj, 1, 2)

    sticks_pairs = pairs_in(vobj.representations["sticks"].indexes)
    assert (1, 2) not in sticks_pairs, "sticks deveria refletir a ligacao removida"
    assert vobj.representations["sticks"].active is True, "active deveria ser preservado"
    assert vobj.representations["lines"].active is False, (
        "lines nao deveria ser ativada a forca so' porque foi recriada"
    )


# ---------------------------------------------------------------------
# Dynamic Bonds (edicao POR FRAME da representacao, nao da topologia)
# ---------------------------------------------------------------------

def test_resolve_frame_arg_forms(make_vismol_object):
    vobj = make_vismol_object(["C", "C"], [0, 1], n_frames=10)

    assert resolve_frame_arg(vobj, None) is None
    assert resolve_frame_arg(vobj, True) == [0]
    assert resolve_frame_arg(vobj, 5) == [5]
    assert resolve_frame_arg(vobj, "1:5") == [1, 2, 3, 4, 5]
    assert resolve_frame_arg(vobj, "all") == list(range(10))


def test_resolve_frame_arg_follows_current_session_frame(make_vismol_object, vm_session):
    vobj = make_vismol_object(["C", "C"], [0, 1], n_frames=10)
    vm_session.set_frame(7)
    assert resolve_frame_arg(vobj, True) == [7]


def test_set_dynamic_bond_order_only_affects_requested_frames(make_vismol_object):
    vobj = make_vismol_object(["C", "C", "N", "H"], [0, 1], n_frames=5)
    vobj.dynamic_bonds = [[0, 1] for _ in range(5)]
    vobj.dynamic_bond_orders = [None] * 5

    n_created = set_dynamic_bond_order(vobj, 2, 3, bond_order=2, frames=[1, 2, 3])

    assert n_created == 3
    assert pairs_in(vobj.dynamic_bonds[0]) == {(0, 1)}
    assert pairs_in(vobj.dynamic_bonds[4]) == {(0, 1)}
    for f in (1, 2, 3):
        assert (2, 3) in pairs_in(vobj.dynamic_bonds[f])
        order = vobj.get_dynamic_bond_order_for_frame(f)
        pair_to_order = dict(zip(
            [tuple(sorted(vobj.dynamic_bonds[f][k:k + 2]))
             for k in range(0, len(vobj.dynamic_bonds[f]), 2)],
            order.tolist(),
        ))
        assert pair_to_order[(2, 3)] == 2, "override manual de ordem deveria ter sido aplicado"


def test_unset_dynamic_bond_only_affects_requested_frame(make_vismol_object):
    vobj = make_vismol_object(["C", "C", "N", "H"], [0, 1], n_frames=5)
    vobj.dynamic_bonds = [[0, 1] for _ in range(5)]
    vobj.dynamic_bond_orders = [None] * 5
    set_dynamic_bond_order(vobj, 2, 3, bond_order=2, frames=[1, 2, 3])

    n_removed = unset_dynamic_bond(vobj, 2, 3, frames=[2])

    assert n_removed == 1
    assert (2, 3) not in pairs_in(vobj.dynamic_bonds[2])
    assert (2, 3) in pairs_in(vobj.dynamic_bonds[1]), "frame 1 nao deveria ter sido afetado"
    assert (2, 3) in pairs_in(vobj.dynamic_bonds[3]), "frame 3 nao deveria ter sido afetado"


def test_dynamic_bonds_do_not_touch_static_topology(make_vismol_object):
    """ Editar Dynamic Bonds e' puramente representacao -- nao pode
    mexer em vismol_object.bonds (a topologia estatica). """
    vobj = make_vismol_object(_BENZENE_SYMBOLS, _BENZENE_BONDS, n_frames=3)
    vobj.dynamic_bonds = [[0, 1] for _ in range(3)]
    vobj.dynamic_bond_orders = [None] * 3
    static_bonds_before = dict(vobj.bonds)

    set_dynamic_bond_order(vobj, 2, 4, bond_order=3, frames=[0, 1, 2])

    assert set(vobj.bonds.keys()) == set(static_bonds_before.keys()), (
        "a topologia estatica (vismol_object.bonds) nao deveria mudar"
    )
