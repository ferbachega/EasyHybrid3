# Modos Normais e Simetria Molecular no EasyHybrid3

## Um guia sobre como o pDynamo3 atribui simetria aos modos vibracionais e por que isso determina a atividade no infravermelho

---

## 1. Objetivo deste tutorial

O EasyHybrid3 permite calcular os modos normais de vibração de um sistema (via `Simulate → Normal Modes`) e, desde a versão atual, também identificar a **simetria** de cada modo (coluna "Symmetry" na janela de análise `Analysis → Normal Modes`). Este documento explica:

1. O que são modos normais e como o pDynamo3 os calcula;
2. O que é um grupo pontual de simetria e como o pDynamo3 o determina automaticamente a partir da geometria 3D;
3. Como cada modo vibracional recebe um rótulo de simetria (uma representação irredutível, como `ag`, `bu`, `a`, `e`, etc.);
4. **Por que a simetria de um modo determina se ele é ativo no infravermelho (IR)** — a regra de seleção espectroscópica — com números reais calculados pelo próprio pDynamo3 para confirmar a teoria.

Todo o texto é escrito do ponto de vista de **como o pDynamo3 implementa isso**, não de uma exposição abstrata de teoria de grupos. Os exemplos usam o *n*-butano (RM1/MNDO, o mesmo sistema usado nos testes internos do EasyHybrid3) nas suas duas conformações — *anti* e *gauche* — porque elas têm grupos pontuais diferentes (C2h e C2, respectivamente) e ilustram bem o contraste entre "com" e "sem" centro de inversão.

---

## 2. Modos normais: o que o pDynamo3 calcula

Na aproximação harmônica, a energia potencial perto de um mínimo (ou ponto estacionário) é aproximada por uma forma quadrática nas coordenadas cartesianas de deslocamento. O pDynamo3 (`pSimulation.NormalModes.NormalModes_SystemGeometry`) faz o seguinte:

1. Calcula a Hessiana (segunda derivada da energia) — analítica quando o modelo de energia oferece, ou numérica por diferenças finitas do gradiente (`ObjectiveFunction.NumericalHessian`) caso contrário;
2. Faz o **mass-weighting** da Hessiana: `H'_ij = H_ij / sqrt(m_i * m_j)`, transformando o problema para coordenadas que já incorporam a massa de cada átomo;
3 Remove (ou projeta para fora) os 6 graus de liberdade de translação e rotação de corpo rígido (`RemoveRotationTranslation`), que não são vibrações reais;
4. Diagonaliza a Hessiana mass-weighted: os autovalores `λ_i` dão as frequências (`ω_i = sqrt(|λ_i|)`, convertidas para cm⁻¹ pelo fator `_To_Wavenumbers`) e os autovetores dão os **vetores de deslocamento** de cada modo (depois "un-mass-weighted" de volta para coordenadas cartesianas reais).

O resultado fica em `system.scratch.nmState` — um objeto `NormalModeState` com `frequencies` (array de frequências, incluindo eventuais modos com frequência negativa/imaginária, sinal de que a estrutura não é um mínimo real) e `modes` (a matriz de autovetores).

Isso é exatamente o que a janela **Simulate → Normal Modes** do EasyHybrid3 executa; o resultado é salvo como uma trajetória `.ptGeo` por modo (`modeN.ptGeo`), com um arquivo `frequency.log` dentro de cada uma.

---

## 3. Grupos pontuais: por que a geometria importa

Um **grupo pontual de simetria** é o conjunto de operações geométricas (rotações, reflexões, inversão, rotações impróprias) que, aplicadas à molécula, produzem uma estrutura indistinguível da original — mantendo pelo menos um ponto fixo no espaço (por isso "pontual": nada se translada).

O pDynamo3 **não pede que você informe o grupo pontual** — ele o determina automaticamente a partir da geometria 3D e dos números atômicos, via `pScientific.Symmetry.Find3DGraphPointGroup`. O algoritmo (em `PointGroupFinder.py`):

1. Translada a estrutura para o centro de massa (ponderado pelas massas atômicas);
2. Classifica os átomos em grupos simetricamente equivalentes (`PostulateSymmetryEquivalentNodes`);
3. Usa o tensor de momento de inércia para restringir as possibilidades (molécula linear, esférica, simétrica, assimétrica);
4. Procura eixos de rotação de ordem superior e operações de ordem 2 (C2, inversão *i*, planos de reflexão σ) que efetivamente mapeiam a estrutura nela mesma, dentro de uma tolerância numérica.

Isso é feito **sobre a estrutura vibracionalmente otimizada** (a mesma para a qual os modos normais foram calculados) — então o grupo pontual detectado reflete a conformação exata que você tem em mãos, não uma classificação abstrata da "molécula" em geral.

### Exemplo real: as duas conformações do butano

| Conformação | Grupo pontual detectado pelo pDynamo3 | Elementos de simetria |
|---|---|---|
| *anti* (trans, cadeia estendida) | **C2h** | E, C2, i (centro de inversão), σh |
| *gauche* | **C2** | E, C2 (só o eixo de rotação, sem centro de inversão) |

Essa diferença — ter ou não um centro de inversão — é o ponto central deste tutorial, como veremos na Seção 6.

---

## 4. Tabelas de caracteres: de onde vêm os rótulos `ag`, `bu`, `a`, `b`...

Cada grupo pontual tem um número finito de **representações irredutíveis** (IRs) — os "tipos de simetria" possíveis para qualquer propriedade da molécula (um modo vibracional, um orbital molecular, um estado eletrônico). Cada IR tem um **caráter** (um número) associado a cada operação de simetria do grupo, organizados na **tabela de caracteres**.

No pDynamo3, essas tabelas **não estão embutidas no código** — são arquivos de dados YAML em `$PDYNAMO3_PARAMETERS/pointGroups/*.yaml`, um por grupo pontual. Por exemplo, o arquivo `C2h.yaml` (usado para a conformação *anti* do butano):

```yaml
Label : C2h
Symmetry Operations :
    C2    :  1
    E     :  1
    i     :  1
    sigma :  1
Character Symmetry Operations : [ E, C2, i ]
Character Table :
    Ag : [ 1.000000 ,  1.000000 ,  1.000000 ]
    Au : [ 1.000000 ,  1.000000 , -1.000000 ]
    Bg : [ 1.000000 , -1.000000 ,  1.000000 ]
    Bu : [ 1.000000 , -1.000000 , -1.000000 ]
```

E o arquivo `C2.yaml` (conformação *gauche*, sem centro de inversão):

```yaml
Label : C2
Character Symmetry Operations : [ E, C2 ]
Character Table :
    A : [ 1.000000 ,  1.000000 ]
    B : [ 1.000000 , -1.000000 ]
```

Repare que **o grupo C2h tem o dobro de representações do grupo C2** — exatamente porque cada representação de C2 (A, B) se desdobra em duas quando o centro de inversão é adicionado: uma versão que não muda de sinal sob inversão (sufixo **g**, de *gerade*, "par" em alemão) e uma que muda de sinal (sufixo **u**, de *ungerade*, "ímpar"). Isso não é uma coincidência de nomenclatura — é a origem exata da regra de exclusão mútua que veremos na Seção 6.

O pDynamo3 carrega esses arquivos via `PointGroup.FromYAML`, guardando a tabela de caracteres como uma matriz (`Array`) indexada por `[IR, operação de simetria]`.

---

## 5. Como cada modo recebe seu rótulo de simetria

Uma vez que o grupo pontual e sua tabela de caracteres são conhecidos, o pDynamo3 precisa descobrir **a que representação irredutível cada modo vibracional pertence**. A função responsável é `pSimulation.NormalModes.NormalModes_IrreducibleRepresentations`, que o EasyHybrid3 chama automaticamente logo depois de `NormalModes_SystemGeometry` (em `p_methods/normal_modes.py`).

O método, em linhas gerais:

1. Para cada operação de simetria do grupo (por exemplo, E, C2, i), calcula o **caráter observado** de cada modo: aplica a operação de simetria geometricamente ao vetor de deslocamento do modo (usando o mapeamento de átomos da operação, `operation.mapping`) e projeta o resultado de volta no vetor original (`vector1.Dot(vector2)`, em `_CharacterFunction`);
2. Compara esse conjunto de caracteres observados com cada linha da tabela de caracteres do grupo pontual, dentro de uma tolerância (`_CharacterMatchTolerance = 0.1`);
3. Atribui a IR cuja linha bate com os caracteres observados.

Quando dois ou mais modos têm frequências muito próximas (dentro de `degeneracyTolerance`, 5 cm⁻¹ por padrão no EasyHybrid3), o pDynamo3 os trata como um possível conjunto degenerado e tenta encontrar uma **combinação** de IRs que explique os caracteres do bloco inteiro — é por isso que às vezes você vê rótulos como `ag/bu`: dois modos que ficaram acidentalmente muito próximos em frequência sem serem realmente relacionados por simetria (degenerescência acidental, não uma degenerescência exigida pelo grupo).

Um rótulo `?` aparece quando nenhuma IR (nem combinação) explica os caracteres observados dentro da tolerância — geralmente ruído numérico da Hessiana (especialmente comum nos "modos" de frequência próxima de zero, que deveriam ser translação/rotação pura, mas raramente saem perfeitamente limpos de um cálculo numérico).

**No EasyHybrid3, o rótulo final é escrito em minúsculas** (`ag`, não `Ag`) — convenção espectroscópica padrão para diferenciar rótulos de modos vibracionais (minúsculo) de rótulos de estados eletrônicos (maiúsculo, `Ag`), embora a tabela de caracteres do pDynamo3 internamente use maiúsculas.

---

## 6. A regra de seleção: por que a simetria decide quem aparece no infravermelho

### 6.1 O critério

Uma transição vibracional fundamental (0 → 1) só absorve luz infravermelha se o **momento de dipolo de transição** for diferente de zero:

```
⟨0| μ |1⟩ ≠ 0
```

onde `μ = (μx, μy, μz)` é o operador de momento de dipolo. Como o estado fundamental vibracional (v=0) é sempre totalmente simétrico (transforma-se como a representação totalmente simétrica do grupo, `Ag` em C2h ou `A` em C2), essa integral só é diferente de zero se a representação do modo vibracional (`Γ(Q)`) for **igual à representação de pelo menos uma componente do momento de dipolo** (x, y ou z) — porque o produto de duas representações iguais sempre contém a representação totalmente simétrica, e o produto de representações diferentes, em geral, não.

Em outras palavras: **um modo só é ativo no IR se ele se transforma da mesma forma que x, y ou z transformam sob as operações de simetria do grupo pontual.**

O pDynamo3 não precisa que você consulte uma tabela de caracteres externa para saber a que representação x, y e z pertencem em cada grupo — ele resolve isso **numericamente**, calculando a intensidade de IR de verdade por diferenças finitas: `pSimulation.NormalModes.NormalModes_InfraredIntensities`. Para cada átomo livre e cada direção cartesiana, desloca a estrutura por ±δ, recalcula o momento de dipolo do modelo QC (`system.DipoleMoment()`), obtém a derivada do dipolo por diferenças finitas, projeta essa derivada nos vetores de modo (`dipoleDerivativesQ.MatrixMultiply(modes, dipoleDerivatives, xTranspose=True)`) e eleva ao quadrado — literalmente calculando `|∂μ/∂Q|²` para cada modo, escalado para km/mol (`_IntensityScalingFactor = 42.2561`).

Ou seja: a intensidade de IR não é um "efeito colateral" derivado da simetria — é uma quantidade calculada diretamente da física (derivada do dipolo). Mas a simetria **prevê exatamente quais modos vão dar zero**, e isso serve como uma verificação cruzada poderosa do cálculo, além de uma explicação física do "porquê".

*(Nota: `NormalModes_InfraredIntensities` está disponível na biblioteca pDynamo3 mas ainda não está exposta em nenhuma janela do EasyHybrid3 — hoje você só vê os rótulos de simetria, não a intensidade numérica. Veja a Seção 8.)*

### 6.2 Caso especial: moléculas com centro de inversão — a regra de exclusão mútua

Isso é onde a distinção `g`/`u` da Seção 4 se torna concreta. Sob a operação de inversão *i*, as coordenadas cartesianas (x, y, z) trocam de sinal — ou seja, **x, y e z sempre se transformam como uma representação `u` (ungerade)** em qualquer grupo com centro de inversão. Consequência direta:

> **Em uma molécula centrossimétrica, nenhum modo `g` (gerade) pode ser ativo no infravermelho — só modos `u` (ungerade) podem.**

Essa é a **regra de exclusão mútua**: em moléculas com centro de inversão, um modo nunca é simultaneamente ativo no infravermelho e no Raman (a polarizabilidade, que rege a atividade Raman, transforma-se como x², y², z², xy, xz, yz — funções sempre `g`). Modos `g` só aparecem no Raman; modos `u` só aparecem no IV.

### 6.3 Confirmação numérica real (butano *anti*, C2h)

Rodei o cálculo completo (`NormalModes_SystemGeometry` + `NormalModes_IrreducibleRepresentations` + `NormalModes_InfraredIntensities`) no *n*-butano na conformação *anti* (RM1/MNDO) para confirmar a regra com números de verdade, não apenas a teoria:

| Modo | Freq. (cm⁻¹) | Simetria | Intensidade IV (km/mol) |
|---:|---:|:---:|---:|
| 10 | 463.6 | **ag** | **0.0000** |
| 13 | 931.1 | **ag** | **0.0000** |
| 17 | 1155.4 | **ag** | **0.0000** |
| 21 | 1205.3 | **ag** | **0.0000** |
| 8 | 192.9 | **bg** | **0.0000** |
| 12 | 887.2 | **bg** | **0.0000** |
| 16 | 1127.6 | **bg** | **0.0000** |
| 7 | 180.9 | au | 0.0097 |
| 11 | 791.7 | au | **4.4647** |
| 14 | 991.7 | au | 1.3213 |
| 9 | 303.0 | bu | 0.2590 |
| 15 | 1019.1 | bu | **6.4903** |
| 23 | 1291.3 | bu | 4.2063 |

**Todo modo `ag` ou `bg` (gerade) deu intensidade exatamente zero (ou ruído numérico ≤ 0.0005). Todo modo `au` ou `bu` (ungerade) deu intensidade real, mensurável.** Essa separação perfeita, obtida de um cálculo de dipolo totalmente independente do cálculo de simetria, é a confirmação prática da regra de exclusão mútua — exatamente como a teoria de grupos prevê, sem nenhum ajuste manual.

### 6.4 E quando não há centro de inversão? (butano *gauche*, C2)

No grupo C2 (sem *i*), **não existe a distinção g/u** — então a regra de exclusão mútua simplesmente não se aplica, e, em princípio, tanto modos `a` quanto `b` podem ser ativos no IV (a que exatamente x, y, z pertencem depende de como os eixos são definidos e teria que ser verificado com `NormalModes_InfraredIntensities`, já que a tabela de caracteres do pDynamo3 não anota "x/y/z" nem "Rx/Ry/Rz" ao lado de cada IR — apenas os valores de caráter numéricos). Isso é consistente com o que vimos na prática: a coluna de simetria da conformação *gauche* mostra só `a`/`b`, e nenhum deles pode ser descartado como "proibido" por simetria antes de calcular a intensidade real.

---

## 7. Resumo prático: como interpretar a coluna "Symmetry"

1. **Abra `Analysis → Normal Modes`**, selecione o sistema/objeto e clique "Import" apontando para as pastas `modeN.ptGeo` que você quer inspecionar. A coluna "Symmetry" mostra o rótulo de cada modo importado.
2. **Verifique o grupo pontual no log do job** (`Simulate → Normal Modes` grava uma linha `Point group: C2h` no log, via Process Manager) — os rótulos de simetria só fazem sentido junto com o grupo pontual detectado.
3. **Se o grupo tem centro de inversão** (contém "h" no nome, como C2h, D2h, Oh, ou o próprio Ci) e a maioria dos seus modos de interesse aparece com sufixo `g`, é bem provável que boa parte deles seja invisível num espectro de infravermelho real — mas ainda pode aparecer no Raman.
4. **Rótulos combinados (`ag/bu`)** indicam duas frequências muito próximas que a rotina não conseguiu separar como pertencentes a uma única IR — normalmente degenerescência acidental (não exigida pela simetria do grupo), especialmente comum em modos de frequência muito baixa (torções, modos quase livres).
5. **Rótulos `?`** geralmente aparecem nos modos de frequência próxima de zero (resíduos de translação/rotação não perfeitamente removidos) — não costuma indicar um problema real com os modos de frequência mais alta.

---

## 8. Limitações e próximos passos

- A **aproximação harmônica** é a base de tudo isso: os rótulos de simetria e a própria existência de modos normais bem definidos dependem de a superfície de energia potencial ser bem aproximada por uma forma quadrática ao redor do ponto calculado. Isso é menos confiável para modos de frequência muito baixa (torções de baixa barreira, por exemplo).
- A **detecção do grupo pontual é numérica e depende de tolerância geométrica**: pequenas distorções (de uma otimização não totalmente convergente, por exemplo) podem levar o pDynamo3 a detectar um grupo pontual de simetria mais baixa do que a simetria "ideal" da molécula (ex.: C1 em vez de C2v). Vale sempre conferir a linha "Point group" no log antes de confiar nos rótulos.
- **Intensidades de infravermelho reais** (`NormalModes_InfraredIntensities`) já existem no pDynamo3 e foram usadas para validar este tutorial, mas **ainda não estão expostas em nenhuma janela do EasyHybrid3** — hoje é preciso rodá-las via script Python separado, como fizemos aqui. Uma extensão natural da janela de análise seria adicionar essa coluna também, usando exatamente a mesma função.
- O mesmo mecanismo de simetria (`Find3DGraphPointGroup` + `IdentifyIrreducibleRepresentations`) também é usado pelo pDynamo3 para identificar a simetria de **orbitais moleculares** e de **estados eletrônicos** (`OrbitalSymmetries`/`StateSymmetries`, em `examples/pSimulation/DetermineSymmetries.py`) — fora do escopo deste tutorial, mas usa exatamente a mesma infraestrutura.

---

## 9. Referências

- Wilson, E. B.; Decius, J. C.; Cross, P. C. *Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra*. McGraw-Hill, 1955.
- Cotton, F. A. *Chemical Applications of Group Theory*, 3rd ed. Wiley, 1990.
- Código-fonte consultado diretamente: `pDynamo3/pSimulation/NormalModes.py`, `pDynamo3/pScientific/Symmetry/PointGroup.py`, `pDynamo3/pScientific/Symmetry/PointGroupFinder.py`, `pDynamo3/parameters/pointGroups/*.yaml`.
