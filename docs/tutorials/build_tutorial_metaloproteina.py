#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Gera docs/tutorials/tutorial_metaloproteina_dyff.html (arquivo unico, imagens embutidas)
#  a partir das capturas em img_metaloproteina/.
#
#      python3 docs/tutorials/build_tutorial_metaloproteina.py
#

import base64
import os

HERE = os.path.dirname ( os.path.abspath ( __file__ ) )
IMG = os.path.join ( HERE, "img_metaloproteina" )
OUT = os.path.join ( HERE, "tutorial_metaloproteina_dyff.html" )


def img ( name, caption, width = "100%" ):
    data = base64.b64encode ( open ( os.path.join ( IMG, name ), "rb" ).read ( ) ).decode ( )
    return ( '<figure><img src="data:image/png;base64,{}" style="max-width:{}" alt="{}">'
             '<figcaption>{}</figcaption></figure>' ).format ( data, width, caption.split ( "." )[0], caption )


CSS = """
:root{--bg:#fff;--fg:#1d2329;--muted:#5b6670;--line:#e2e6ea;--card:#f6f8fa;--accent:#2b6cb0;--warn:#c05621;--ok:#2f855a}
@media (prefers-color-scheme:dark){:root{--bg:#15191d;--fg:#e6e9ec;--muted:#9aa5b1;--line:#2c3339;--card:#1d2328;--accent:#7fb2ea;--warn:#f6ad55;--ok:#68d391}}
body{background:var(--bg);color:var(--fg);font:16px/1.65 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;margin:0}
main{max-width:980px;margin:0 auto;padding:28px 16px 80px}
h1{font-size:1.8em;line-height:1.25;margin:.2em 0 .4em}h2{margin-top:2.4em;border-bottom:2px solid var(--line);padding-bottom:.25em}
h3{margin-top:1.6em}.muted{color:var(--muted)}
figure{margin:1.2em 0;text-align:center}figure img{border:1px solid var(--line);border-radius:6px;box-shadow:0 2px 10px rgba(0,0,0,.15)}
figcaption{font-size:.9em;color:var(--muted);margin-top:.4em;text-align:left}
.box{background:var(--card);border-left:5px solid var(--accent);padding:12px 16px;border-radius:6px;margin:1.2em 0}
.box.warn{border-left-color:var(--warn)}.box.ok{border-left-color:var(--ok)}
.box h4{margin:.1em 0 .4em}
.step{display:inline-block;background:var(--accent);color:#fff;border-radius:50%;width:1.7em;height:1.7em;text-align:center;line-height:1.7em;margin-right:.4em;font-weight:600}
table{border-collapse:collapse;width:100%;font-size:.93em;margin:.8em 0;display:block;overflow-x:auto}
th,td{border:1px solid var(--line);padding:6px 9px;text-align:left;vertical-align:top}th{background:var(--card)}
code{background:var(--card);padding:1px 5px;border-radius:4px;font-size:.9em}
pre{background:var(--card);padding:12px;border-radius:6px;overflow-x:auto;font-size:.85em;line-height:1.45}
.menu{font-family:system-ui;background:var(--card);border:1px solid var(--line);border-radius:4px;padding:1px 6px;white-space:nowrap}
nav{background:var(--card);border-radius:6px;padding:10px 16px;margin:1em 0}nav ol{margin:.3em 0;padding-left:1.4em}
"""

html = """<!doctype html><html lang="pt-BR"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Tutorial: Metaloproteína com DYFF</title><style>{css}</style></head><body><main>
<p class="muted">EasyHybrid 3 · tutorial · preparo de sistemas com DYFF refinado</p>
<h1>Preparando uma metaloproteína para simulação com DYFF: a anidrase carbônica II</h1>
<p>Este tutorial mostra, passo a passo, como transformar uma estrutura cristalográfica de uma metaloenzima em um sistema pronto para
simulação no EasyHybrid, usando o campo de força <b>DYFF refinado</b>: biblioteca de resíduos ajustada a dados DFT, e o sítio metálico
parametrizado na hora com cálculos xTB. Além do <i>como</i>, explicamos o <i>porquê</i> de cada decisão. Em metaloproteínas, a maioria dos
erros de preparo não vem do programa, mas de escolhas químicas: estado de oxidação, spin, protonação dos ligantes e modelo do metal.</p>

<nav><b>Roteiro</b><ol>
<li><a href="#sistema">O sistema de estudo</a></li>
<li><a href="#fundamentos">Fundamentos: por que metais são difíceis</a></li>
<li><a href="#p1">Importar a estrutura</a></li>
<li><a href="#p2">Abrir o preparo DYFF</a></li>
<li><a href="#p3">Entender e ajustar a janela de preparo</a></li>
<li><a href="#protonacao">Fundamentos: protonação dos ligantes do metal</a></li>
<li><a href="#outros">Outros casos: capas covalentes, cobre, resíduos truncados</a></li>
<li><a href="#p4">Executar: o que acontece por dentro</a></li>
<li><a href="#p5">Ler o resultado</a></li>
<li><a href="#p6">Inspecionar o sítio do zinco</a></li>
<li><a href="#p7">Validar: o sítio se mantém?</a></li>
<li><a href="#limites">Cuidados e limitações</a></li>
<li><a href="#script">Fazendo o mesmo por script</a></li>
<li><a href="#refs">Referências</a></li></ol></nav>

<h2 id="sistema">1. O sistema de estudo</h2>
<p>A <b>anidrase carbônica II humana</b> (CA II) catalisa a hidratação do CO<sub>2</sub> (CO<sub>2</sub> + H<sub>2</sub>O ⇌ HCO<sub>3</sub><sup>−</sup> + H<sup>+</sup>),
uma das reações mais rápidas da bioquímica. O sítio ativo tem um íon <b>Zn<sup>2+</sup></b> coordenado por três histidinas (His94, His96 e His119)
e por uma molécula de água ou hidróxido, que é o nucleófilo da reação. Usamos a estrutura <b>PDB 3KS3</b>, resolvida a 0,9 Å, que inclui
480 águas cristalográficas e uma molécula de glicerol (GOL) usada como crioprotetor.</p>
<table><tr><th>Item</th><th>Valor</th></tr>
<tr><td>Arquivo</td><td><code>tests/data/metalloproteins/3KS3.pdb</code> (a base de testes inclui outras 13 proteínas: <code>README.md</code> e <code>manifest.json</code> na mesma pasta)</td></tr>
<tr><td>Sítio</td><td>Zn<sup>2+</sup> · His94 NE2 2,03 Å · His96 NE2 2,04 Å · His119 ND1 2,02 Å · H<sub>2</sub>O 1,90 Å</td></tr>
<tr><td>Pré-requisitos</td><td>EasyHybrid com o <code>xtb</code> configurado (versão 6.7.1 ou posterior para usar g-xTB) e RDKit no Python do EasyHybrid</td></tr>
<tr><td>Tempo</td><td>cerca de 1 min 10 s para o preparo completo com GFN2-xTB (≈ 4 min com g-xTB)</td></tr></table>

<h2 id="fundamentos">2. Fundamentos: por que metais são difíceis para campos de força</h2>
<p>Um campo de força clássico descreve a molécula com termos simples: ligações e ângulos como molas, torções como séries de cossenos, e
cargas fixas com Lennard-Jones para as interações não ligadas. Isso funciona muito bem para a química orgânica, que se repete: um grupo
metila é parecido em qualquer molécula. Metais de transição quebram essas premissas:</p>
<ul>
<li><b>a ligação metal–ligante é parcialmente covalente</b> e transfere carga. O Zn "+2" de um sítio real carrega bem menos que +2;</li>
<li><b>a geometria depende do metal, do estado de oxidação e do spin</b>. O Zn<sup>2+</sup> (d<sup>10</sup>) não tem preferência eletrônica e
costuma ser tetraédrico; o Cu<sup>2+</sup> (d<sup>9</sup>) sofre distorção de Jahn–Teller; o Fe<sup>3+</sup> de alto spin tem 5 elétrons
desemparelhados;</li>
<li><b>os parâmetros não são transferíveis</b>: um Zn–N(His) num sítio catalítico não é igual ao de um sítio estrutural.</li></ul>
<p>Há três estratégias clássicas de modelagem:</p>
<table><tr><th>Modelo</th><th>Como funciona</th><th>Vantagens</th><th>Problemas</th></tr>
<tr><td><b>Não ligado (iônico)</b></td><td>O metal é uma carga pontual com Lennard-Jones; nenhuma ligação com a proteína.</td><td>Simples; os ligantes podem trocar.</td>
<td>A geometria de coordenação não é garantida: ligantes escapam, e o número de coordenação muda.</td></tr>
<tr><td><b>Ligado (<i>bonded</i>)</b></td><td>Ligações, ângulos e cargas metal–ligante explícitos, parametrizados por cálculo quântico do sítio (a ideia do MCPB.py do AMBER).</td>
<td>Mantém a geometria do sítio; cargas mais realistas.</td><td>Não permite troca de ligantes; precisa de cálculo por sítio.</td></tr>
<tr><td><b>QM/MM</b></td><td>O sítio é tratado por mecânica quântica, e o resto pelo campo de força.</td><td>O mais fiel; permite reações.</td><td>Muito mais caro.</td></tr></table>
<div class="box"><h4>O que o EasyHybrid faz</h4>
O preparo DYFF usa o <b>modelo ligado</b> como padrão para metais de transição e o <b>iônico</b> para Mg<sup>2+</sup> e Ca<sup>2+</sup>, que
são íons "duros" cuja interação com oxigênios é essencialmente eletrostática. O sítio ligado é parametrizado em três partes:
<ol><li>um <b>subsistema com capas</b> é recortado em volta do metal: metal, resíduos doadores e vizinhança de duas ligações, com valências
fechadas por H;</li>
<li>a <b>Hessiana</b> desse subsistema é calculada com GFN2-xTB ou g-xTB, e as constantes de ligação e ângulo saem pelo método de
<b>Seminario modificado</b>, que projeta a matriz de segundas derivadas nas direções de cada ligação e ângulo;</li>
<li>as <b>cargas</b> são CM5 (GFN1-xTB) do mesmo subsistema.</li></ol>
Os demais resíduos vêm da <b>biblioteca de resíduos</b>, cujos ângulos e torções foram ajustados às forças DFT de cerca de 330 dipeptídeos. Na avaliação
com outros dipeptídeos, essa biblioteca ficou mais precisa nas forças que o OPLS-AA.</div>

<h2 id="p1"><span class="step">1</span>Importar a estrutura</h2>
<p>Abra o PDB por <span class="menu">File › Import Data...</span> (ou pelo ícone de abrir). A figura mostra a CA II recém-importada:
cadeia em linhas, águas cristalográficas como cruzes vermelhas.</p>
{img01}
<div class="box warn"><h4>Por que ainda não dá para simular</h4>
Uma estrutura cristalográfica <b>não tem hidrogênios</b>. Os raios X praticamente não os enxergam, salvo em resolução sub-atômica. O
arquivo também não traz <b>ordens de ligação</b> nem, em geral, registros CONECT. O EasyHybrid desenha ligações pela distância, só para
visualização, mas o sistema pDynamo não tem conectividade química. Sem H, sem ordens de ligação e sem estados de protonação, nenhum campo
de força pode ser atribuído. É isso que o preparo resolve.</div>

<h2 id="p2"><span class="step">2</span>Abrir o preparo DYFF</h2>
<p>Use <span class="menu">Extras › Prepare DYFF System...</span>. O mesmo preparo também abre sozinho se você escolher
<span class="menu">System › Add › Force Field › DYFF + residue library...</span> num sistema que ainda não tem ligações.</p>
{img02}
<div class="box"><h4>Por que não o "Prepare OPLS System"?</h4>
O preparo OPLS do EasyHybrid só conhece os íons do OPLS-AA (Na<sup>+</sup>, K<sup>+</sup>, Cl<sup>−</sup>, Mg<sup>2+</sup>, Ca<sup>2+</sup>).
Com Zn, Cu ou Fe ele para com "<i>no OPLS ion parameters</i>". O preparo DYFF reaproveita a análise da estrutura (resíduos, pontes
dissulfeto, terminais) e a construção a partir da biblioteca de componentes do PDB, mas não tipa nada com OPLS. Por isso aceita metais.</div>

<h2 id="p3"><span class="step">3</span>Entender e ajustar a janela de preparo</h2>
{img03}
<p>Cada parte da janela, de cima para baixo:</p>
<h3>Resumo da estrutura</h3>
<p><i>"740 residues (257 amino acid, 2 unsupported, 481 water)"</i>. Os 2 resíduos "sem componente" são o Zn e o glicerol; nenhum dos
dois existe na biblioteca de componentes do pDynamo, e os dois recebem tratamento próprio. Não há quebras de cadeia nem pontes
dissulfeto (a CA II não tem cistinas).</p>
<h3>Sítios metálicos</h3>
<p>Uma linha por metal. A coluna <b>Donors</b> lista os átomos doadores encontrados e suas distâncias: <i>HOH.402 O 1,90; HIS.119 ND1 2,02;
HIS.94 NE2 2,03; HIS.96 NE2 2,04</i>. É exatamente o sítio descrito na literatura. A detecção usa a soma dos raios covalentes mais 0,5 Å
(para Zn–N, até 2,43 Å), o que inclui ligações longas e exclui contatos de segunda esfera. As colunas <b>Model</b>, <b>Charge</b> e
<b>Multiplicity</b> são editáveis com clique duplo.</p>
<div class="box"><h4>Fundamento: carga e multiplicidade</h4>
A <b>carga</b> é o estado de oxidação formal do metal. A <b>multiplicidade</b> 2S+1 vem do número de elétrons desemparelhados: o xTB
recebe <code>--uhf</code> = multiplicidade − 1. Errar a multiplicidade muda a estrutura eletrônica do cálculo e, com ela, as constantes
e as cargas do sítio.
<table><tr><th>Íon</th><th>Configuração</th><th>Multiplicidade típica em proteínas</th></tr>
<tr><td>Zn<sup>2+</sup>, Cd<sup>2+</sup></td><td>d<sup>10</sup></td><td>1</td></tr>
<tr><td>Cu<sup>2+</sup></td><td>d<sup>9</sup></td><td>2</td></tr>
<tr><td>Cu<sup>+</sup></td><td>d<sup>10</sup></td><td>1</td></tr>
<tr><td>Fe<sup>3+</sup> alto spin (rubredoxina oxidada)</td><td>d<sup>5</sup></td><td>6</td></tr>
<tr><td>Fe<sup>2+</sup> alto spin</td><td>d<sup>6</sup></td><td>5</td></tr>
<tr><td>Ni<sup>2+</sup> octaédrico</td><td>d<sup>8</sup></td><td>3</td></tr>
<tr><td>Mn<sup>2+</sup></td><td>d<sup>5</sup></td><td>6</td></tr></table>
Os valores padrão da janela seguem essa tabela. Confira sempre com a literatura do seu sistema: estados de baixo spin existem, por
exemplo em hemes.</div>
<h3>Ligantes</h3>
<p>O glicerol recebe hidrogênios e ordens de ligação pelo Open Babel, no pH indicado. Ele é tratado como o sítio metálico: cargas CM5,
Hessiana e torções. Glicerol é um crioprotetor e não pertence à proteína em solução; para simulações de produção, considere removê-lo
antes, no próprio arquivo PDB.</p>
<h3>pH, água implícita, método, escala CM5</h3>
<ul>
<li><b>pH</b>: usado só para os hidrogênios dos ligantes. A protonação dos resíduos da proteína vem da análise da estrutura e, nos
ligantes do metal, do próprio sítio (próxima seção).</li>
<li><b>Implicit water (ALPB)</b>: os cálculos de cargas e Hessiana do subsistema são feitos em água implícita. Um sítio carregado em fase
gasosa tem uma estrutura eletrônica irrealista: carboxilatos, por exemplo, ficam instáveis. Mantenha marcado.</li>
<li><b>Hessian method</b>:
<ul><li><i>GFN2-xTB</i> é rápido e aceita água implícita;</li>
<li><i>g-xTB</i> é mais preciso, mas roda pelo tblite, que não tem solvente implícito. Para sítios carregados, prefira GFN2-xTB;
para o cobre, o g-xTB descreve melhor a ligação axial fraca com a Met (seção 7);</li>
<li><i>DFT (ORCA)</i> já aparece, mas ainda não está disponível.</li></ul></li>
<li><b>CM5 charge scale</b>: fator que multiplica as cargas mantendo a carga total de cada resíduo. Contra dados DFT em fase gasosa, valores
entre 0,7 e 0,8 deram as melhores energias conformacionais. Isso favorece cargas menores do que o ideal em solução, onde ainda não há
avaliação. Mantenha 1,0, a menos que você teste.</li></ul>

<h2 id="protonacao">4. Fundamentos: protonação dos ligantes do metal</h2>
<p>Este é o erro mais comum, e o mais silencioso, no preparo de metaloproteínas. Um átomo que doa um par de elétrons ao metal
<b>não pode estar protonado</b>:</p>
<ul>
<li><b>Histidina</b>: o anel imidazol tem dois nitrogênios, ND1 (δ) e NE2 (ε). Em pH neutro, só um deles carrega H. Se o Zn se liga pelo
NE2, o H tem de estar no ND1 (tautômero <b>HID</b>, "Delta Protonated"), e vice-versa (<b>HIE</b>). Um preparo que escolhe o tautômero
"mais comum" sem olhar o metal coloca um H apontando para o zinco;</li>
<li><b>Cisteína</b>: ligada ao metal, é <b>tiolato</b> (S<sup>−</sup>, sem H), como nos sítios Zn–Cys<sub>4</sub> e Fe–Cys<sub>4</sub>;</li>
<li><b>Asp/Glu</b>: carboxilatos desprotonados.</li></ul>
<p>O preparo DYFF aplica essas regras automaticamente aos resíduos dos sítios <i>ligados</i>. Na CA II, a mensagem final informa:
<i>"HIS.119 Epsilon Protonated, HIS.94 Delta Protonated, HIS.96 Delta Protonated"</i>. A His119 coordena pelo ND1, por isso o H vai para
o NE2; a His94 e a His96 coordenam pelo NE2, por isso o H vai para o ND1.</p>
<div class="box warn"><h4>A água do zinco: água ou hidróxido?</h4>
Na CA II, a água ligada ao Zn tem pK<sub>a</sub> ≈ 7: em pH fisiológico, coexistem Zn–OH<sub>2</sub> e Zn–OH<sup>−</sup>, e é o hidróxido que
ataca o CO<sub>2</sub>. O preparo mantém a molécula como <b>água</b>, sem ligação ao Zn: o DYFF não tem tipo atômico para um oxigênio de
água com três ligações. A interação Zn···O fica puramente eletrostática, e a mensagem final avisa isso. Se o seu estudo depende do
estado catalítico (Zn–OH<sup>−</sup>), o hidróxido precisa ser montado à mão (por exemplo, editando o sistema no Builder): nesta versão,
o preparo não converte a água coordenada em hidróxido nem a liga ao metal.</div>

<h2 id="outros">5. Outros casos que a janela trata</h2>
<h3>Capas covalentes: dedo de zinco (PDB 1ZNF)</h3>
{img04}
<p>Estruturas de RMN de peptídeos costumam ter <b>capas</b> ACE (acetil) e NH2 (amida) ligadas covalentemente à cadeia. Esses resíduos não
estão na biblioteca de componentes, e o preparo não sabe ligá-los à proteína. Tratados como moléculas soltas, ficariam a 1,3 Å da cadeia
sem ligação, uma sobreposição catastrófica. A janela os detecta e oferece removê-los: a cadeia fica com terminais livres (NH<sub>3</sub><sup>+</sup>
e COO<sup>−</sup>). O sítio, Zn com 2 Cys e 2 His, é o clássico sítio estrutural dos dedos de zinco: as Cys serão tiolatos.</p>
<h3>Cobre de camada aberta: plastocianina (PDB 1PLC)</h3>
{img05}
<p>O cobre "azul" (tipo 1) tem duas His, uma Cys e uma <b>Met axial longa</b> (Cu–SD 2,82 Å). A janela já propõe multiplicidade 2 para
Cu<sup>2+</sup>. Atenção: o modelo ligado trata a interação longa com a Met como uma ligação comum e a encurta na minimização (2,82 → 2,38 Å
com GFN2-xTB; 2,55 Å com g-xTB). Para cobre, prefira <b>g-xTB</b> e verifique essa distância.</p>
<h3>Resíduos truncados</h3>
{img08}
<p>Cadeias laterais desordenadas costumam ser omitidas no depósito: Lys e Glu de superfície são os casos típicos. O preparo pode
<b>completar</b> esses resíduos. Os átomos que faltam são colocados geometricamente, com comprimentos e ângulos padrão (tetraédricos
nos carbonos sp<sup>3</sup>) e conformação estendida. É uma aproximação razoável para começar, não uma predição de rotâmero: minimize e,
se o resíduo for importante para o seu estudo, investigue os rotâmeros. A alternativa clássica, truncar para Ala ou Gly, muda a química.</p>

<h2 id="p4"><span class="step">4</span>Executar: o que acontece por dentro</h2>
<p>Clique em <b>Prepare and apply DYFF</b>. Na CA II, o processo leva cerca de 70 s com GFN2-xTB e a barra de status acompanha cada etapa:</p>
<table><tr><th>#</th><th>Etapa</th><th>Detalhe</th></tr>
<tr><td>1</td><td>Construção</td><td>Hidrogênios e ordens de ligação a partir da biblioteca de componentes do PDB, com as protonações escolhidas. As
ligações Zn–N(His) entram na conectividade.</td></tr>
<tr><td>2</td><td>Reparo de hidrogênios</td><td>O construtor do pDynamo coloca um H quase em cima do átomo pai quando os vizinhos são
quase coplanares, o que acontece em geometrias cristalográficas distorcidas. Esses H são refeitos com geometria ideal, e águas que se
chocam são reorientadas.</td></tr>
<tr><td>3</td><td>Tipo do metal</td><td>O DYFF tipa átomos por padrões, e muitos metais não têm padrão para a coordenação encontrada. O
preparo atribui o tipo pelo elemento, pela carga e pelo número de coordenação: aqui, <code>Zn:Tet:+2</code>.</td></tr>
<tr><td>4</td><td>Biblioteca de resíduos</td><td>254 de 259 resíduos reconhecidos pelo grafo de ligações. As três His do sítio ficam
<b>pendentes</b> de propósito: uma ligação externa ao metal impede o casamento com a biblioteca, e esses resíduos são parametrizados
junto com o Zn. O glicerol também fica fora.</td></tr>
<tr><td>5</td><td>Refinamento dos grupos fora da biblioteca</td><td>Dois grupos: (i) Zn + 3 His, num subsistema de 83 átomos com capas;
(ii) glicerol. Em cada um: cargas CM5 (GFN1-xTB, ALPB), Hessiana GFN2-xTB e constantes de Seminario <b>calibradas</b> ao DFT pelo ambiente
químico (78 termos). As 17 ligações giráveis recebem torções da biblioteca DFT.</td></tr>
<tr><td>6</td><td>Não ligados</td><td>H polares ganham um Lennard-Jones pequeno, pois com cargas parciais ficariam dentro da parede
repulsiva dos aceptores. As 481 águas recebem cargas TIP3P e o LJ do oxigênio da TIP3P.</td></tr>
<tr><td>7</td><td>Modelo</td><td>O DYFF é atribuído com todas as edições, e o modelo não ligado é definido.</td></tr></table>

<h2 id="p5"><span class="step">5</span>Ler o resultado</h2>
{img06}
<p>O sistema novo, <b>3KS3_DYFF</b>, aparece na árvore de objetos. O original fica preservado, aqui ocultado. A mensagem completa da
barra de status, na aba <i>Status</i> do painel inferior, é a seguinte:</p>
<pre>MM model 'DYFF (dyff-1.0)' assigned to '3KS3_DYFF'. (7 guanidinium group(s) corrected ...)
(5512 partial charges from the object: Builder or residue library; 1412 polar hydrogens with a small
Lennard-Jones radius; 481 water oxygens with the TIP3P Lennard-Jones.)
Residue library: 254/259 residues; total charge +1.000.
Outside the library (GFN2-xTB): 67 bonds, 110 angles from the Hessian; 17 of 17 rotatable bonds from
the DFT torsion library.
Prepared: protonation set by the metal sites: HIS.119 Epsilon Protonated, HIS.94 Delta Protonated,
HIS.96 Delta Protonated; 1 bonded metal site(s) refined with GFN2-xTB.
WARNING: A:ZN.262: the coordinated water A:HOH.402 is kept without a metal-oxygen bond ...</pre>
<table><tr><th>Trecho</th><th>Significado</th><th>O que conferir</th></tr>
<tr><td><i>guanidinium group(s) corrected</i></td><td>Correção automática da tipagem do DYFF nas argininas.</td><td>—</td></tr>
<tr><td><i>5512 partial charges</i></td><td>Todos os átomos têm carga.</td><td>Deve ser igual ao número de átomos.</td></tr>
<tr><td><i>total charge +1.000</i></td><td>Soma das cargas formais da proteína com Zn<sup>2+</sup>.</td><td>Um valor inteiro e coerente com a
protonação. Não neutralizado: se for simular em solvente explícito, adicione contraíons.</td></tr>
<tr><td><i>254/259 residues</i></td><td>Os 5 restantes são as 3 His do sítio, a água ligada e o glicerol.</td><td>Muitos resíduos não
reconhecidos indicam nomes ou protonações incomuns.</td></tr>
<tr><td><i>17 of 17 rotatable bonds</i></td><td>Todas as torções do sítio e do glicerol vieram da biblioteca DFT.</td><td>Ligações sem
correspondência ficam com as regras do DYFF.</td></tr>
<tr><td><i>WARNING ... coordinated water</i></td><td>A água do Zn interage só por eletrostática (seção 4).</td><td>Decida entre água e
hidróxido.</td></tr></table>

<h2 id="p6"><span class="step">6</span>Inspecionar o sítio do zinco</h2>
{img07}
<p>Os resíduos do sítio e a água estão em <i>sticks</i>, e o Zn é o pequeno átomo cinza no centro. Observe:</p>
<ul>
<li>as <b>três ligações Zn–N</b>, que fazem parte do modelo MM;</li>
<li>o <b>H de cada His</b> está no nitrogênio que <i>não</i> coordena: o átomo ligado ao Zn não tem H;</li>
<li>a <b>água</b> aparece junto ao Zn. A visualização desenha contatos curtos, mas no modelo MM não há ligação Zn–O.</li></ul>
<p>Para reproduzir esta visão: selecione os resíduos do sítio, mostre-os como <i>sticks</i> pelo menu de representação da seleção
(botão direito na área gráfica) e centralize a visão na seleção.</p>

<h2 id="p7"><span class="step">7</span>Validar: o sítio se mantém?</h2>
<p>A validação mais direta é <b>minimizar</b> o sistema e conferir se as distâncias metal–doador continuam próximas das experimentais:</p>
<ol>
<li><span class="menu">Simulate › Geometry optimization</span>: algumas centenas a milhares de passos;</li>
<li><span class="menu">Analysis › Distances, Angles and Dihedrals</span>, ou a seleção por clique (pk1/pk2), para medir Zn–NE2(His94),
Zn–NE2(His96), Zn–ND1(His119) e Zn–O(água).</li></ol>
<p>Fizemos esse teste para os três modelos: ligado com GFN2-xTB, ligado com g-xTB e iônico. Foram 1000 passos de LBFGS, no vácuo, sem as
águas que não se ligam ao metal:</p>
{img09}
<table><tr><th>Doador</th><th>Cristal</th><th>Ligado, GFN2-xTB</th><th>Ligado, g-xTB</th><th>Iônico</th></tr>
<tr><td>His94 NE2</td><td>2,03</td><td>1,94</td><td>1,93</td><td>2,21</td></tr>
<tr><td>His96 NE2</td><td>2,04</td><td>1,91</td><td>1,95</td><td>2,12</td></tr>
<tr><td>His119 ND1</td><td>2,02</td><td>1,98</td><td>1,96</td><td><b>3,68</b></td></tr>
<tr><td>H<sub>2</sub>O O</td><td>1,90</td><td>2,40</td><td>2,42</td><td>2,21</td></tr></table>
<div class="box ok"><h4>Interpretação</h4>
No <b>modelo ligado</b>, as três His ficam a no máximo 0,13 Å das distâncias do cristal, com os dois métodos. No <b>modelo iônico</b>, a
His119 se solta do zinco (3,68 Å): o sítio se desfaz. A água, sem ligação no modelo ligado, se afasta para 2,4 Å; essa é a consequência
da escolha discutida na seção 4. O mesmo padrão aparece em outros sítios da base de testes: na rubredoxina (Fe<sup>3+</sup>–Cys<sub>4</sub>),
o maior desvio é 0,03 Å no modelo ligado e 0,14 Å no iônico.</div>

<h2 id="limites">8. Cuidados e limitações</h2>
<ul>
<li><b>Sem validação contra DFT para metais.</b> As bibliotecas DFT, de torções e de calibração, vêm de moléculas orgânicas do QCArchive;
os termos em torno do metal vêm só do xTB. Para conclusões quantitativas sobre o sítio, compare com QM/MM ou DFT.</li>
<li><b>Ligações fracas</b>, como Cu–S(Met) axial, ficam fortes demais no modelo ligado.</li>
<li><b>A água coordenada</b> não é ligada (seção 4).</li>
<li><b>Sítios completados por simetria cristalográfica</b>, como o Zn da insulina 4INS, têm parte dos ligantes em cópias de simetria
que não estão no arquivo. Com menos de dois doadores, o sítio vira iônico, com aviso. Monte o conjunto biológico antes.</li>
<li><b>Cofatores com metal</b>, como heme e clusters Fe–S, ainda não são suportados pelo preparo DYFF.</li>
<li><b>Água explícita</b>: as cargas TIP3P sobre termos ligados do DYFF não formam um modelo de água validado. Equilibre e confira a
densidade antes de produzir.</li>
<li><b>Modelo ligado não permite troca de ligantes.</b> Para estudar a entrada de substrato no sítio ou reações, use QM/MM.</li></ul>

<h2 id="script">9. Fazendo o mesmo por script</h2>
<p>O preparo também pode ser usado diretamente em Python, o que é útil para preparar muitas estruturas:</p>
<pre>import pBabel
from pBabel import ImportSystem
from pdynamo import dyff_prep as DP, residue_library as RL

system = ImportSystem("3KS3.pdb")
info = DP.analyze(system)                      # residues, metal sites, ligands, truncated, covalent
DP.remove_residues(info, info["covalent"])     # caps / covalent ligands (none in 3KS3)
for site in info["metals"]:                    # check/edit: model, charge, multiplicity
    print(site["label"], site["element"], site["model"], site["charge"], site["multiplicity"],
          [(info["analysis"].residues[j].label, name, round(d, 2)) for j, name, d in site["donors"]])
DP.apply_site_protonation(info)                # His tautomers, thiolates, carboxylates
prepared, report = DP.build(info, label="3KS3_DYFF")
plan = RL.plan(prepared)                       # residue library matching
rest = RL.rest_parameters(prepared, plan["rest_atoms"], solvent="water", method="gfn2",
                          site_states=DP.site_groups(prepared, report))
# then: charges + edits on the object and define_MMModel('DYFF') -- what the dialog does
# (main_window._run_dyff_library), with report["type_overrides"] for the metal types.</pre>
<p>Para validar uma série de proteínas como na seção 7: <code>python3 fraglib/fraglib/run_metal_sites.py --ids 3KS3 1IRO --methods gfn2 gxtb ionic</code>.</p>

<h2 id="refs">10. Referências</h2>
<ol>
<li>J. M. Seminario, <i>Int. J. Quantum Chem.</i> 60, 1271 (1996); A. E. A. Allen, M. C. Payne, D. J. Cole, <i>J. Chem. Theory Comput.</i> 14, 274 (2018) (Seminario modificado).</li>
<li>C. Bannwarth, S. Ehlert, S. Grimme, <i>J. Chem. Theory Comput.</i> 15, 1652 (2019) (GFN2-xTB).</li>
<li>S. Ehlert et al., <i>J. Chem. Theory Comput.</i> 17, 4250 (2021) (ALPB).</li>
<li>A. V. Marenich et al., <i>J. Chem. Theory Comput.</i> 8, 527 (2012) (cargas CM5).</li>
<li>P. Li, K. M. Merz Jr., <i>J. Chem. Inf. Model.</i> 56, 599 (2016) (MCPB.py: o modelo ligado para metaloproteínas).</li>
<li>D. N. Silverman, R. McKenna, <i>Acc. Chem. Res.</i> 40, 669 (2007) (mecanismo da anidrase carbônica).</li>
<li>Material Suplementar do EasyHybrid: <code>docs/dyff_refinement_SI.html</code> (validação das bibliotecas DYFF);
fluxograma de preparo: <code>docs/dyff_preparo_sistema.html</code>.</li></ol>
</main></body></html>"""

html = html.replace ( "{css}", CSS )
for entry in (
    ( "img01", "01_imported.png", "Figura 1. A CA II (3KS3) recém-importada. A árvore de objetos mostra o sistema <i>0 - 3KS3</i>; a barra inferior informa "
                                  "2537 átomos (sem hidrogênios) e <i>Force Field: None</i>.", "100%" ),
    ( "img02", "02_menu_extras.png", "Figura 2. O menu <i>Extras</i> com o item <i>Prepare DYFF System...</i>.", "40%" ),
    ( "img03", "03_dialog.png", "Figura 3. A janela de preparo para a CA II: um sítio de Zn com quatro doadores, modelo ligado, carga 2, multiplicidade 1; "
                                "glicerol como ligante; GFN2-xTB em água implícita.", "100%" ),
    ( "img04", "04_dialog_zinc_finger.png", "Figura 4. Dedo de zinco 1ZNF: capas ACE e NH2 ligadas covalentemente (marcadas para remoção) e o sítio Zn–Cys<sub>2</sub>His<sub>2</sub>.", "100%" ),
    ( "img05", "05_dialog_plastocyanin.png", "Figura 5. Plastocianina 1PLC: Cu<sup>2+</sup> com His<sub>2</sub>, Cys e a Met axial a 2,82 Å; multiplicidade 2 (d<sup>9</sup>).", "100%" ),
    ( "img08", "08_dialog_truncated.png", "Figura 6. Uma ubiquitina com as cadeias laterais de duas Lys removidas: a janela lista os átomos que faltam e oferece completá-los, com aviso.", "100%" ),
    ( "img06", "06_result.png", "Figura 7. Resultado: o sistema <i>1 - 3KS3_DYFF</i>, em cartoon, e a mensagem completa na aba Status.", "100%" ),
    ( "img07", "07_site.png", "Figura 8. O sítio do zinco no sistema preparado: His94, His96 e His119 ligadas ao Zn, com os H nos nitrogênios que não coordenam; a água coordenada em vermelho e branco.", "100%" ),
    ( "img09", "09_site_distances.png", "Figura 9. Distâncias Zn–doador no cristal e após minimização, nos três modelos.", "100%" ) ):
    key, name, caption, width = entry
    html = html.replace ( "{" + key + "}", img ( name, caption, width ) )

with open ( OUT, "w" ) as handle: handle.write ( html )
print ( OUT, "{:.1f} MB".format ( os.path.getsize ( OUT ) / 1e6 ) )
