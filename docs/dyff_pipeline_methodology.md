# Pipeline de montagem de sistemas no DYFF-1.0eh

*Resumo dos princípios adotados na preparação de sistemas DYFF-1.0eh (Seminário modificado, cargas CM5, solvente implícito, biblioteca de resíduos). 2026-11-16.*

## 1. Visão geral (dado bruto → sistema pronto)

```
PDB / MOL2 / coordenadas
  → percepção de ligações (quando não há bonds no arquivo)
  → Prepare DYFF System (hidrogênios, protonação, sítios metálicos)
  → biblioteca de resíduos (match por isomorfismo de grafo)
  → átomos fora da biblioteca → xTB (CM5 + Seminário + torções)
  → DefineMMModel (DYFF + edits por termo)
```

- **Percepção de ligação:** `pdynamo.bond_order_perception` (Wang & Case) para PDB sem conectividade; MOL2 já traz ordens de ligação (`mol2_import.py`, com fallback via obabel quando a leitura nativa falha).
- **Prepare DYFF System** (`dyff_prep.py`): constrói hidrogênios, decide estado de protonação, detecta sítios metálicos por distância (soma dos raios covalentes + 0.5 Å), remove capas covalentes (ACE/NH2), completa resíduos truncados.
- **Biblioteca de resíduos** (`residue_library.py`): compara cada resíduo contra ~26 templates por isomorfismo de grafo. Resíduo que bate = usa parâmetros **pré-computados por DFT** (cargas CM5, bonds/angles de Seminário, torções). Isso cobre os 20 aminoácidos padrão e variantes comuns — é o caminho barato e validado.
- **Átomos fora da biblioteca** (`rest_atoms`): ligantes, cofatores, resíduos não-padrão, sítios metálicos, Cys em dissulfeto. Esses vão para xTB sob demanda, em subsistemas capeados com H (`subsystem.py`), processados em paralelo por grupo conectado.

## 2. Princípio de base do DYFF

DYFF não tem arquivo de parâmetros — é "UFF-like": tipo de átomo vem só de **elemento + hibridização + ordem de ligação**, e os termos (bond/angle/dihedral/out-of-plane) são gerados por **regras geométricas genéricas**. O que fazemos é **corrigir esses termos genéricos** com dados de química quântica, termo a termo, guardando os valores como *edits* por par/trio de átomos (`dyff_edits.py`, `MMModelDYFFEdited`) em vez de reescrever o campo de força inteiro. Cada edit é amarrado a uma "assinatura de topologia" — se a conectividade local mudar, o edit é descartado em vez de aplicado errado.

## 3. Cargas parciais — CM5

- **Padrão:** CM5 (Charge Model 5) sobre **GFN1-xTB**, com **solvente implícito ALPB (água)** por padrão (`solvent="water"` em `rest_charges`/`xtb_cm5`). É o mesmo esquema usado na biblioteca de resíduos (validado contra PDBCharges, MAE 0.010 e em 1ALF).
- **Alternativa para ligantes carregados:** RESP sobre GFN2-xTB ou HF/6-31G* (`ligand_charges="resp-xtb"/"resp-hf"`), escolhido quando CM5 não é suficiente (zwitteríons, grupos muito polarizados).
- **Escala opcional:** `cm5_scale` reescala dentro de cada grupo/resíduo mantendo o total (`q' = q_mean + s(q - q_mean)`); testado em dipeptídeos SPICE (gas-phase, 0.7–0.8 ótimo), **não validado em fase condensada** — por isso o padrão é 1.0 e só mexemos nisso quando testado.

## 4. Bonds e angles — Seminário modificado

- Vêm do **Hessiano de xTB** (`seminario.py`), não de tabelas: para cada par/trio de átomos "fora da biblioteca", o bloco da matriz Hessiana projetado na direção da ligação/ângulo vira constante de força (método de Seminário modificado, Allen et al. 2018 — o mesmo princípio do MCPB.py do AMBER para metais).
- Métodos disponíveis: **GFN2-xTB** (padrão, aceita solvente implícito) ou **g-xTB** (mais caro, mais preciso, só fase gasosa — tblite não tem solvente implícito). Disponível tanto no "Prepare DYFF System" quanto no refinamento via janela DYFF Parameters.
- **Guard-rails:** `MIN_BOND_K`/`MIN_ANGLE_K` — constante não-física (ex.: carboxilato em fase gasosa dá k negativo) faz o termo manter o valor padrão do DYFF em vez do valor do xTB; fica registrado em `skipped`.
- **Calibração opcional:** quando existe biblioteca de Hessianos DFT por ambiente químico (`bond_library.py`), os valores de xTB são corrigidos (`shift`, `ratio`) antes do guard-rail.
- **Correção de 2026-11-16 (parte 1 — âncoras):** o subsistema capeado roda `xtb --ohess` (otimiza geometria **e depois** tira o Hessiano). Sem nada fixo, um sítio multidentado (metal com 4 doadores) relaxa livremente no vácuo e perde a geometria real da proteína. Fixamos só os átomos de fronteira (ligados a uma capa de H). Essa parte é real e reproduzível, mas **não bastou**: o sítio continuou distorcendo depois da otimização do sistema completo.
- **Correção de 2026-11-16 (parte 2 — ângulo doador-metal-doador é descartado do ajuste):** o DYFF já tem um mecanismo próprio para o ângulo doador-M-doador (geometria ideal percebida da estrutura + constante de força estilo UFF pelas cargas efetivas, calibrado contra tmQM — `dyff_edits.py`'s `_preset_metal_geometry`/`_metal_angles`). Um edit por termo do Hessian para o MESMO trio sobrepõe esse mecanismo, e o ajuste do Hessian do fragmento capeado é consistentemente mole demais para esse ângulo específico quando o resto da proteína volta a entrar em jogo (1ZNF: k ajustado ~95-150 kJ/mol/rad², sítio caiu para 80-134° depois de otimizar o sistema completo, em 3 execuções independentes). Agora `seminario.drop_metal_centered_angles` descarta esse edit (mantém ligação e os ângulos centrados no doador, ex. Cys CB-SG-Zn, que não competem com nada). Revalidado pelo pipeline real (1ZNF): 113.8°/105.1° (Cys-Zn-Cys/His-Zn-His) depois de otimizar, idêntico ao cristal. Ver `fraglib/reports/fraglib_status.md` seção 17 e 17.1.

## 5. Sítios metálicos

- Detectados por distância, não por nome de resíduo. Dois modelos:
  - **bonded**: ligações explícitas metal–doador, tipado via `MMModelDYFFEdited.typeOverrides`, refinado com o Hessiano (Seminário) junto dos resíduos doadores.
  - **ionic**: íon livre (sem ligação explícita) quando doadores insuficientes — mais simples, mas perde ligantes fracos (ex.: Met axial em plastocianina tipo 1).
- 17 metais têm estado de oxidação/multiplicidade/raio padrão (`METALS`); metal fora dessa lista agora é **sinalizado** (`unrecognized_metals`) em vez de virar íon silencioso sem aviso.

## 6. Torções

- Preferência por **biblioteca DFT** (QCArchive TorsionDrive, B3LYP-D3BJ/DZVP), buscada por ambiente químico (exato → 2 ligações → 1 ligação ao redor da ligação) — tão bom quanto ou melhor que escanear a própria molécula.
- Ligação flexível sem cobertura na biblioteca: **scan relaxado em xTB** (24 pontos) + ajuste de série de cosseno.
- Torção só é usada **junto** com bonds/angles de Seminário — a biblioteca foi validada em cima dessa base, não isolada.

## 7. Out-of-plane / impropers

Resultado negativo documentado: tentamos extrair a rigidez de piramidalização do próprio Hessiano, mas não deu ganho — os impropers continuam pelas regras genéricas do DYFF, não por xTB.

## 8. Não-ligados

- LJ estilo OPLS (11 tipos), não os parâmetros nativos do DYFF.
- Hidrogênios polares têm LJ suavizado (`soften_polar_hydrogens`) — evita over-binding em pontes de H.
- Água tem LJ calibrado à parte (constantes próprias, não genéricas).

## 9. Validação (resumo)

- 928 Hessianos DFT (QCArchive): Seminário+xTB reduz RMSE de frequência de ~180 para ~59–126 cm⁻¹ dependendo do nível.
- SPICE (dipeptídeos capeados): forças/energias contra DFT.
- Metaloproteínas (14 PDBs): desvio de distância metal-doador ≤ 0.1–0.5 Å na maioria dos sítios — mas essa tabela **não** media ângulo, que é onde o bug de 2026-11-16 apareceu.
- Proteína-ligante: PLA15, PL-REX, OpenFF Protein-Ligand Benchmark, PLINDER (não PDBbind/Binding MOAD).

## 10. Pontos em aberto

- Cofatores metálicos complexos (heme, clusters Fe-S) não suportados.
- DFT (ORCA) como alternativa ao xTB está na interface mas **desligado** (não implementado).
- A correção de âncoras do item 4 só foi revalidada no 1ZNF — os outros 5 sítios metálicos da tabela de validação e os fragmentos de torção da biblioteca de resíduos passam pelo mesmo mecanismo e ainda não foram re-rodados.
