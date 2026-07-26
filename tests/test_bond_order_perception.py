"""
Testes para a percepcao de ordem de ligacao por casamento maximo
(vismol.core.bond_order_perception.perceive_bond_order_for_pairs_pure).

Cobre os bugs reais encontrados e corrigidos durante o desenvolvimento
(ver o historico de commits em src/graphics_engine/src/vismol/core/
bond_order_perception.py para o contexto completo de cada um):

  1) GABEDIT_MAX_VALENCE tinha um bug de maiusculas/minusculas que fazia
     a valencia maxima de qualquer elemento de 2 letras (Cl, Na, Fe, ...)
     ser ignorada silenciosamente.
  2) A atribuicao de duplas era feita por um guloso de uma passada so',
     que dependia da ORDEM em que os pares apareciam no array de entrada
     -- podendo dar uma estrutura quimicamente invalida (ex.: 1,3-
     butadieno com uma extremidade sem a dupla que precisava) dependendo
     so' de como o arquivo/parser ordenou os bonds.
  3) Atomos de fronteira QC/MM (Dynamic Bonds) tem o grau local
     subestimado, porque a ligacao para a regiao MM nao aparece no
     subconjunto de pares -- sem correcao, isso promove indevidamente
     uma ligacao vizinha a dupla.

Esses testes nao precisam do VismolObject completo (nem das fixtures em
conftest.py) -- perceive_bond_order_for_pairs_pure e' uma funcao PURA,
sem estado, exatamente para ser testavel isoladamente.
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

from vismol.core.bond_order_perception import (  # noqa: E402
    perceive_bond_order_for_pairs_pure,
    GABEDIT_MAX_VALENCE,
)


def _orders_by_pair(flat_pairs, order):
    n_bonds = len(flat_pairs) // 2
    return {
        tuple(sorted((flat_pairs[2 * k], flat_pairs[2 * k + 1]))): order[k]
        for k in range(n_bonds)
    }


def _assert_valid_valences(symbols, flat_pairs, order):
    """ Nenhum atomo pode ter a soma das ordens de ligacao MAIOR que sua
    valencia maxima -- viola-la e' o tipo de erro que a versao antiga
    (guloso de uma passada) podia produzir silenciosamente. """
    n_bonds = len(flat_pairs) // 2
    valence = {}
    for k in range(n_bonds):
        i, j = flat_pairs[2 * k], flat_pairs[2 * k + 1]
        valence[i] = valence.get(i, 0) + order[k]
        valence[j] = valence.get(j, 0) + order[k]
    for atom, v in valence.items():
        mx = GABEDIT_MAX_VALENCE.get(symbols[atom], 4)
        assert v <= mx, (
            f"atomo {atom} ({symbols[atom]}) com valencia {v}, "
            f"acima do maximo permitido ({mx})"
        )


# ---------------------------------------------------------------------
# 1) Bug de case do GABEDIT_MAX_VALENCE (Cl, Na, Fe, ...)
# ---------------------------------------------------------------------

def test_chlorine_never_promoted_to_double():
    # Cl deveria ter valencia maxima 1 -- nunca promovida a dupla, mesmo
    # que o carbono vizinho tenha folga de valencia (Cl-CH3).
    symbols = ["Cl", "C", "H", "H", "H"]
    bonds = [0, 1, 1, 2, 1, 3, 1, 4]
    order = perceive_bond_order_for_pairs_pure(symbols, bonds)
    pair_order = _orders_by_pair(bonds, order)
    assert pair_order[(0, 1)] == 1
    _assert_valid_valences(symbols, bonds, order)


def test_sodium_chloride_stays_single():
    symbols = ["Na", "Cl"]
    bonds = [0, 1]
    order = perceive_bond_order_for_pairs_pure(symbols, bonds)
    assert order[0] == 1


# ---------------------------------------------------------------------
# 2) Sistema saturado (etano) -- nao deve inventar duplas
# ---------------------------------------------------------------------

def test_ethane_all_single_bonds():
    symbols = ["C", "C", "H", "H", "H", "H", "H", "H"]
    bonds = [0, 1, 0, 2, 0, 3, 0, 4, 1, 5, 1, 6, 1, 7]
    order = perceive_bond_order_for_pairs_pure(symbols, bonds)
    assert all(o == 1 for o in order)


# ---------------------------------------------------------------------
# 3) Butadieno: order-independent (bug do guloso de uma passada)
# ---------------------------------------------------------------------

def test_butadiene_order_independent():
    symbols = ["C", "C", "C", "C", "H", "H", "H", "H", "H", "H"]
    bonds_natural_order = [0, 1, 1, 2, 2, 3, 0, 4, 0, 5, 1, 6, 2, 7, 3, 8, 3, 9]
    # Ligacao central aparece PRIMEIRO no array -- e' a ordem que
    # quebrava o algoritmo guloso antigo.
    bonds_middle_first = [1, 2, 0, 1, 2, 3, 0, 4, 0, 5, 1, 6, 2, 7, 3, 8, 3, 9]

    order_a = perceive_bond_order_for_pairs_pure(symbols, bonds_natural_order)
    order_b = perceive_bond_order_for_pairs_pure(symbols, bonds_middle_first)
    pa = _orders_by_pair(bonds_natural_order, order_a)
    pb = _orders_by_pair(bonds_middle_first, order_b)

    assert pa[(0, 1)] == 2
    assert pa[(1, 2)] == 1
    assert pa[(2, 3)] == 2
    assert pb == pa, "resultado nao deveria depender da ordem dos bonds de entrada"

    _assert_valid_valences(symbols, bonds_natural_order, order_a)
    _assert_valid_valences(symbols, bonds_middle_first, order_b)


# ---------------------------------------------------------------------
# 4) Benzeno: alternancia correta, tambem order-independent
# ---------------------------------------------------------------------

def test_benzene_alternates_correctly():
    symbols = ["C"] * 6 + ["H"] * 6
    bonds_1 = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0,
               0, 6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11]
    bonds_2 = [3, 4, 0, 1, 1, 2, 2, 3, 4, 5, 5, 0,
               0, 6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11]

    order_1 = perceive_bond_order_for_pairs_pure(symbols, bonds_1)
    order_2 = perceive_bond_order_for_pairs_pure(symbols, bonds_2)
    p1 = _orders_by_pair(bonds_1, order_1)
    p2 = _orders_by_pair(bonds_2, order_2)

    n_double_1 = sum(1 for k, v in p1.items() if v == 2 and k[0] < 6 and k[1] < 6)
    n_double_2 = sum(1 for k, v in p2.items() if v == 2 and k[0] < 6 and k[1] < 6)
    assert n_double_1 == 3
    assert n_double_2 == 3

    _assert_valid_valences(symbols, bonds_1, order_1)
    _assert_valid_valences(symbols, bonds_2, order_2)


# ---------------------------------------------------------------------
# 5) Naftaleno (aneis fundidos): valido e deterministico
# ---------------------------------------------------------------------

def test_naphthalene_fused_rings_deterministic():
    symbols = ["C"] * 10 + ["H"] * 8
    bonds_1 = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 4, 6, 6, 7, 7, 8, 8, 9, 9, 5,
               0, 10, 1, 11, 2, 12, 3, 13, 6, 14, 7, 15, 8, 16, 9, 17]
    bonds_2 = [0, 10, 0, 1, 1, 11, 1, 2, 2, 12, 2, 3, 3, 13, 3, 4, 4, 6, 4, 5,
               6, 14, 6, 7, 7, 15, 7, 8, 8, 16, 8, 9, 9, 17, 9, 5, 5, 0]

    order_1 = perceive_bond_order_for_pairs_pure(symbols, bonds_1)
    order_2 = perceive_bond_order_for_pairs_pure(symbols, bonds_2)
    p1 = _orders_by_pair(bonds_1, order_1)
    p2 = _orders_by_pair(bonds_2, order_2)

    assert p1 == p2, "mesma molecula deveria dar a mesma estrutura de Kekule"
    _assert_valid_valences(symbols, bonds_1, order_1)
    _assert_valid_valences(symbols, bonds_2, order_2)


# ---------------------------------------------------------------------
# 6) Imidazol (anel impar, heteroatomo -- His/purinas)
# ---------------------------------------------------------------------

def test_imidazole_odd_ring_with_heteroatom():
    symbols = ["N", "C", "N", "C", "C", "H", "H", "H"]
    bonds = [0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 1, 5, 3, 6, 4, 7]
    order = perceive_bond_order_for_pairs_pure(symbols, bonds)
    _assert_valid_valences(symbols, bonds, order)

    p = _orders_by_pair(bonds, order)
    n_double = sum(1 for v in p.values() if v == 2)
    assert n_double == 2, "imidazol (5 ligacoes no anel) deve ter exatamente 2 duplas"


# ---------------------------------------------------------------------
# 7) Performance: sistema tipo coroneno (varios aneis fundidos)
# ---------------------------------------------------------------------

def test_coronene_like_system_is_fast():
    import time

    symbols = ["C"] * 24
    ring_center = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    outer = []
    for k in range(6):
        a, b = k, (k + 1) % 6
        o1, o2, o3 = 6 + 3 * k, 6 + 3 * k + 1, 6 + 3 * k + 2
        outer += [(a, o1), (o1, o2), (o2, o3), (o3, b)]

    bonds = []
    for i, j in ring_center + outer:
        bonds += [i, j]

    degree = {}
    for i in range(0, len(bonds), 2):
        degree[bonds[i]] = degree.get(bonds[i], 0) + 1
        degree[bonds[i + 1]] = degree.get(bonds[i + 1], 0) + 1

    hidx = 24
    while len(symbols) <= hidx:
        symbols.append("H")
    for atom in range(24):
        while degree.get(atom, 0) < 3:
            bonds += [atom, hidx]
            while len(symbols) <= hidx:
                symbols.append("H")
            degree[atom] = degree.get(atom, 0) + 1
            hidx += 1

    t0 = time.time()
    order = perceive_bond_order_for_pairs_pure(symbols, bonds)
    dt = time.time() - t0

    assert dt < 2.0, f"deveria terminar rapido (componentes pequenos); levou {dt:.3f}s"
    _assert_valid_valences(symbols, bonds, order)


# ---------------------------------------------------------------------
# 8) Fronteira QC/MM: atomo com ligacao "invisivel" para a regiao MM
# ---------------------------------------------------------------------

def test_qcmm_boundary_bug_reproduction():
    """ Reproduz o bug: sem extra_degree, um N de fronteira (ja saturado
    na valencia real -- 1 ligacao para C, 1 para H, 1 "invisivel" para a
    regiao MM) tem seu grau local subestimado e a ligacao N-C e'
    promovida indevidamente a dupla. """
    symbols = ["N", "C", "H"]
    bonds = [0, 1, 0, 2]  # falta, de proposito, a ligacao N-(regiao MM)
    order = perceive_bond_order_for_pairs_pure(symbols, bonds)
    assert order[0] == 2, "reproducao do bug esperada sem a correcao de fronteira"


def test_qcmm_boundary_fix_with_extra_degree():
    symbols = ["N", "C", "H"]
    bonds = [0, 1, 0, 2]
    order = perceive_bond_order_for_pairs_pure(symbols, bonds, extra_degree={0: 1})
    assert order[0] == 1, "com a correcao, N-C deve ficar simples (N ja saturado)"
    assert order[1] == 1


def test_extra_degree_zero_is_a_no_op():
    symbols = ["N", "C", "H"]
    bonds = [0, 1, 0, 2]
    order_sem = perceive_bond_order_for_pairs_pure(symbols, bonds)
    order_com_zero = perceive_bond_order_for_pairs_pure(symbols, bonds, extra_degree={0: 0})
    assert list(order_sem) == list(order_com_zero)
