# Testes do EasyHybrid3

## Como rodar

```bash
pip install -r requirements-dev.txt --user     # so' adiciona pytest
cd src/graphics_engine && bash install.sh && cd ../..   # compila as extensoes Cython (mesmo comando do install.py)
python3 -m pytest tests/ -v
```

Não é preciso instalar PyGObject/pycairo/GTK3 para rodar esta suíte —
ver "O que esta coberto" abaixo.

## O que está coberto

Esta suíte testa a **lógica de dados** de EasyHybrid/VisMol que não
depende de GTK, OpenGL, ou de um display real:

- **`test_bond_order_perception.py`** — a percepção de ordem de ligação
  por casamento máximo (`vismol.core.bond_order_perception.
  perceive_bond_order_for_pairs_pure`): o bug de case do
  `GABEDIT_MAX_VALENCE`, a dependência de ordem do algoritmo guloso antigo
  (butadieno, benzeno, naftaleno, imidazol), e a correção de fronteira
  QC/MM (`extra_degree`).

- **`test_atom_ops_bonds.py`** — as funções de edição de ligação usadas
  pelo comando de terminal `bond`/`unbond`
  (`gui/windows/builder/atom_ops.py`): `set_bond_order`/`unset_bond`
  (garantindo que editar uma ligação numa estrutura carregada de arquivo
  não apaga as outras — esse foi um bug real, encontrado testando o
  terminal manualmente), o refresh da representação `sticks`, e
  `set_dynamic_bond_order`/`unset_dynamic_bond`/`resolve_frame_arg`
  (edição da representação de Dynamic Bonds por frame).

Isso é possível porque `vismol.core.vismol_object` (e os módulos que
dependem dele) só precisam, para essa parte da lógica, de `numpy` +
`Cython` (compilado, ver abaixo) + `freetype-py` + `PyOpenGL` — **não**
importam `gi`/GTK diretamente. `tests/conftest.py` fornece um
`vm_session`/`vm_config` mínimo (sem nenhuma sessão gráfica real) só com
os atributos que essas funções efetivamente leem.

## O que NÃO está coberto (limitação conhecida)

- Qualquer coisa que dependa de renderização OpenGL de verdade
  (`vismol.libgl.representations`, o próprio desenho na tela) ou da
  interface GTK (`src/gui/windows/setup/easyhybrid_terminal.py` e o
  restante de `src/gui/`) não é exercitada por estes testes — precisaria
  de um display (real ou virtual, ex. Xvfb) e de uma sessão GTK/OpenGL
  completa, o que foge do escopo deste primeiro CI "básico".
- pDynamo3 em si não é uma dependência instalável via pip (ver
  `install.py`/`check_pdynamo()`), então nada que dependa de um sistema
  pDynamo real (ex. `sync_pdynamo_system`) é testado aqui — as funções
  testadas chamam essas integrações por trás de fixtures/dublês que não
  as exercitam de verdade.

Se/quando isso for expandido, os candidatos naturais são: um teste de
sintaxe (`python -m compileall`) para todo o `src/` (já incluído no CI,
ver `.github/workflows/ci.yml`), e, futuramente, um job separado com
Xvfb + GTK3 typelib instalados via `apt-get` para exercitar a parte
gráfica.
