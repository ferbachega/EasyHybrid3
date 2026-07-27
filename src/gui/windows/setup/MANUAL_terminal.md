# Manual do Terminal do EasyHybrid

O terminal do EasyHybrid permite controlar a sessão de visualização por comandos de texto, sem passar pelos menus. Os mesmos comandos também funcionam como API Python para automação.

## Sintaxe

Todo comando segue o formato:

```
comando arg1=valor1 arg2=valor2
```

Os argumentos são nomeados (`chave=valor`), podem vir em qualquer ordem, e os tipos são convertidos automaticamente: `true`/`false` viram booleano, números viram int ou float, o resto é texto. Caminhos com espaço precisam de aspas: `load file="/minha pasta/sis.pdb"`.

Comandos desconhecidos retornam uma mensagem de erro em vez de quebrar a janela, e um argumento errado mostra a forma de uso esperada.

## Teclas úteis

As setas para cima e para baixo navegam pelo histórico de comandos. A tecla Tab completa o nome do comando que você está digitando (se houver só uma opção, completa direto; se houver várias, completa até o prefixo comum). Pressionar Tab duas vezes seguidas lista todas as opções que combinam com o que foi digitado.

## Referência de comandos

### help

Lista todos os comandos disponíveis com uma breve descrição de cada um.

```
help
```

### list

Lista os objetos moleculares carregados, com seus índices. O índice é o número que você usa no argumento `obj=` dos outros comandos.

```
list
```

Saída de exemplo:

```
  [0] proteina.pdb
  [1] ligante.mol2
```

### show e hide

Mostram ou ocultam uma representação. Sem filtros, agem sobre a seleção ativa. Com `obj=` (e opcionalmente os filtros finos), miram apenas aqueles átomos sem alterar a seleção ativa — a ação é pontual.

O argumento `rep=` aceita os tipos suportados pela sessão: `lines`, `sticks`, `spheres`, `dash`, entre outros. Quando omitido, assume `lines`.

```
show rep=sticks
show rep=sticks obj=1
show rep=spheres obj=0 chain=A
show rep=sticks obj=0 resn=HIS name=CA
hide rep=lines obj=0 resi=10-25
hide rep=spheres
```

### select e deselect

O `select` define a seleção ativa de forma persistente: os comandos seguintes que não especificarem `obj=` passam a agir sobre ela, até você fazer um novo `select` ou limpar com `deselect`.

```
select obj=0
select obj=0 chain=A
select obj=0 chain=A resn=HIS
select obj=0 resi=10-30
select obj=0 chain=A resn=HIS name=CA
deselect
```

Fluxo típico — selecionar uma vez e aplicar vários comandos:

```
select obj=0 chain=A
show rep=sticks
hide rep=lines
center
```

### Filtros de seleção

Os filtros são aceitos por `show`, `hide`, `select` e `center`, e se combinam com lógica E (cada filtro adicional restringe mais o resultado).

O filtro `chain` compara o nome da cadeia e aceita um valor único ou uma lista separada por vírgula, como `chain=A` ou `chain=A,B`. O filtro `resi` compara o índice do resíduo e aceita valor único, faixa ou combinação: `resi=45`, `resi=10-20`, `resi=1-5,8,12`. O filtro `resn` compara o nome do resíduo, como `resn=HIS`. O filtro `name` compara o nome do átomo, como `name=CA`.

```
select obj=0 chain=A,B
select obj=0 resi=10-20
select obj=0 resi=1-5,8,12
show rep=sticks obj=0 chain=A resn=TRP
```

### frame, next, prev

Controlam a navegação na trajetória. O `frame` sem argumento mostra o frame atual; com `n=` pula para um frame específico. O `next` e o `prev` avançam e retrocedem um frame.

```
frame
frame n=1000
next
prev
```

### center

Centraliza a câmera num alvo, com uma animação de translação. Sem `obj=`, centraliza na seleção ativa (no centroide dos átomos selecionados). Com `obj=` e sem filtros finos, centraliza no centro de massa do objeto. Com filtros, centraliza no centroide dos átomos que passam no filtro.

```
center
center obj=0
center obj=0 chain=A
center obj=0 resn=HIS name=CA
```

### zoom

Aproxima ou afasta a câmera, equivalente a girar o scroll do mouse. O argumento `dir=` aceita `in` (aproxima) ou `out` (afasta), e `steps=` define quantos passos de scroll aplicar de uma vez (padrão 5).

```
zoom dir=in
zoom dir=out steps=10
```

### axes

Mostra ou oculta os eixos de referência.

```
axes show=true
axes show=false
```

### load

Carrega uma molécula de um arquivo. Caminhos com espaço precisam de aspas.

```
load file=/home/usuario/proteina.pdb
load file="/home/usuario/minha pasta/sistema.xyz"
```

## Receitas práticas

Isolar e destacar uma região de interesse, partindo de um sistema com proteína e ligante:

```
list
hide rep=lines obj=0
show rep=sticks obj=0 chain=A resn=HIS
center obj=0 chain=A resn=HIS
zoom dir=in steps=8
```

Inspecionar uma trajetória passo a passo numa região específica:

```
select obj=0 resi=45-50
show rep=sticks
center
frame n=0
next
next
```

## Uso como API em scripts

Fora do terminal, os mesmos comandos são chamáveis como métodos do objeto `cmd`, com argumentos Python. Isso é útil para automação — por exemplo, varrer uma trajetória de 2000 frames.

```python
cmd.show(rep="sticks", obj=1, chain="A")
cmd.select(obj=0, resi="10-30")
cmd.center(obj=0, chain="A")
cmd.zoom(dir="out", steps=3)

# percorrer a trajetória de 100 em 100 frames
for i in range(0, 2000, 100):
    cmd.frame(n=i)
    # aqui você poderia capturar uma imagem, medir uma distância, etc.
```

A diferença é apenas a forma de escrever: na entrada de texto do terminal usa-se a sintaxe `comando arg=valor`; em script Python usa-se `cmd.comando(arg=valor)`. As duas rotas executam exatamente a mesma lógica.

## Notas

Os filtros `resi` comparam contra o índice interno do resíduo. Se a numeração que você conhece (a do arquivo PDB original) não coincidir com esse índice, uma faixa como `resi=10-25` pode selecionar um intervalo diferente do esperado — vale conferir com um resíduo de número conhecido no primeiro uso.

O comando de cor ainda não está disponível e será adicionado em uma versão futura.
