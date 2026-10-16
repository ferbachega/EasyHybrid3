# Ferramentas de análise de dinâmica molecular — candidatas para o EasyHybrid3

Levantamento do que já existe em `src/gui/windows/analysis/` (RMSD, distância/
ângulo/diedro, WHAM, PES, normal modos, refinamento de energia, alinhamento e
reimaging de trajetória) e do que complementaria essas ferramentas, organizado
por categoria e pensando no perfil do EasyHybrid: QM/MM via pDynamo3 + MD
clássica via NAMD + bicamadas lipídicas.

## Estrutural / distribuição de pares

- **g(r) (RDF)** — parcial (por tipo de átomo/par de seleções) e total; base
  para número de coordenação (integral de g(r) até o primeiro mínimo).
- **Perfil de densidade** — massa/número ao longo de um eixo (essencial para
  bicamada: perfil ao longo de z mostrando cabeça polar, cauda, água).
- **Função de distribuição espacial (SDF)** — versão 3D do g(r), útil para
  camadas de solvatação anisotrópicas.

## Dinâmica / dependente do tempo

- **MSD → coeficiente de difusão** — relação de Einstein, precisa de
  reimaging correto (já existe).
- **Função de autocorrelação de velocidade (VACF)** — densidade de estados
  vibracionais, difusão.
- **Tempo de residência / probabilidade de sobrevivência** — ex. quanto tempo
  uma água fica na primeira esfera de coordenação de um íon.

## Flutuação estrutural / conformacional

- **RMSF** — complemento natural do RMSD que já existe (flutuação por
  átomo/resíduo, identifica regiões flexíveis).
- **Raio de giração (Rg)**.
- **PCA / dinâmica essencial** — modos coletivos de movimento a partir da
  trajetória (diferente dos normais modos que já existem, que são
  harmônicos/estáticos).
- **Clustering conformacional** — agrupar frames por RMSD
  (k-means/hierárquico), achar estruturas representativas.
- **Mapa de contatos** — nativo vs. não-nativo, útil para
  dobramento/desnaturação.

## Ligação de hidrogênio e interações

- **Análise de pontes de H** — contagem, tempo de vida, pares
  doador-aceptor, geometria (distância + ângulo).
- **Pontes salinas / pareamento iônico**.

## Específico de bicamada lipídica

(dado que já trabalhamos UFF + bicamada nesta sessão — ver
`docs/dyff_to_namd_methodology.md`)

- **Área por lipídeo**.
- **Espessura da bicamada** (distância P-P ou equivalente).
- **Parâmetro de ordem das caudas (S_CD)** — o clássico para medir
  ordenação/fluidez de cadeia lipídica.
- **Difusão lateral** — MSD restrito ao plano da membrana.
- **Espectro de ondulação/curvatura da membrana**.

## Superfície e geometria

- **SASA** (área de superfície acessível ao solvente) — nota:
  `surface_analysis_window.py` hoje é *renderização* de isosuperfícies QM
  (orbital/densidade/MEP), não SASA.
- **Detecção de cavidades/bolsões** — relevante para sítios ativos/ligantes.

## Energético / termodinâmico

- **Energia de interação entre grupos** (MM-MM, QM-MM) — decomposição por
  par de seleções.
- **Perfis de energia livre além do WHAM** — MBAR, ou leitura de
  metadinâmica se algum dia entrar.

## Específico de QM/MM

(diferencial do EasyHybrid frente a ferramentas MD genéricas)

- **Evolução de carga/spin da região QM ao longo da trajetória**.
- **Acompanhamento de coordenada de reação** (distância de ligação
  formando/quebrando) — parcialmente coberto por distance/angle/dihedral,
  mas útil um modo "trajetória inteira" dedicado.
- **Momento de dipolo ao longo do tempo** — `System.DipoleMoment()` já
  existe no pDynamo3, só falta expor numa janela/série temporal.

## Sugestão de prioridade

- **RDF + densidade** como par inicial natural — compartilham praticamente
  toda a infraestrutura de binning/PBC.
- **RMSF** logo em seguida — quase gratuito já tendo o RMSD pronto.
