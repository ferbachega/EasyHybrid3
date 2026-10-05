#!/usr/bin/env python3
"""
Exemplo de uso do QCModelNWChem.

Mostra tres cenarios:
  (A) script standalone minimo para TESTAR o modulo com um sistema pequeno
      (o jeito mais rapido de validar os pontos [VALIDATE] contra um NWChem real);
  (B) como o modelo seria criado dentro do EasyHybrid (session.py);
  (C) um sistema QC/MM (regiao quantica + ambiente MM via point charges Bq).

Pre-requisitos:
  - NWChem instalado e a variavel de ambiente definida:
        export PDYNAMO3_NWCHEMCOMMAND=/caminho/para/nwchem
  - PDYNAMO3_SCRATCH definido (ou o modulo cai no tempdir do sistema).
  - QCModelNWChem instalado no pDynamo (pMolecule/QCModel/) -- use o instalador
    do EasyHybrid, ou copie o arquivo e registre no __init__.py.
"""

# ---------------------------------------------------------------------------
# (A) SCRIPT STANDALONE MINIMO  -- teste rapido do modulo
# ---------------------------------------------------------------------------
# Este e o caminho recomendado para VALIDAR o parser: um sistema pequeno
# (agua), energia + gradiente, sem MM. Rode e compare a energia/gradiente com
# uma execucao manual do NWChem para confirmar que o parser le certo.

from pBabel                    import ImportSystem
from pMolecule.QCModel         import QCModelNWChem      # <- vem do pDynamo apos instalar
from pScientific.Geometry3     import Coordinates3
from pSimulation               import ConjugateGradientMinimize_SystemGeometry

# 1) carregar (ou construir) um sistema com uma regiao QC definida
system = ImportSystem ( "agua.xyz" )          # qualquer sistema pequeno
system.electronicState.charge       = 0
system.electronicState.multiplicity = 1

# 2) criar o modelo NWChem
qcModel = QCModelNWChem.WithOptions (
    method     = "dft"      ,   # dft | scf | mp2 ...
    functional = "b3lyp"    ,   # xc do bloco dft
    basis      = "6-31G*"   ,   # base (biblioteca do NWChem)
    charge     = 0          ,
    # memory   = "2000 mb"  ,   # opcional
    # keywords = [ "iterations 200", "convergence energy 1e-7" ],  # linhas extra no bloco dft
    scratch    = None       ,   # None -> resolve PDYNAMO3_SCRATCH / tempdir
)

# opcional: apontar o executavel explicitamente (senao usa PDYNAMO3_NWCHEMCOMMAND)
# qcModel.SetCommand ( "/opt/nwchem/bin/nwchem" )

# 3) associar o modelo ao sistema e calcular a energia
system.DefineQCModel ( qcModel )
energy = system.Energy ( )                    # roda o NWChem, le a energia
print ( "Energia (kJ/mol):", energy )

# 4) as propriedades extraidas ficam em target.scratch.nwchemOutputData
data = system.scratch.nwchemOutputData
print ( "Energia (Hartree) :", data.get ( "Energy" ) )
print ( "Convergiu?        :", data.get ( "Is Converged" ) )
print ( "Terminou normal?  :", data.get ( "Normal Termination" ) )
print ( "Cargas Mulliken   :", data.get ( "Mulliken Charges" ) )
print ( "Dipolo            :", data.get ( "Dipole" ) )

# 5) uma otimizacao de geometria (usa o gradiente do NWChem)
ConjugateGradientMinimize_SystemGeometry ( system,
                                           maximumIterations    = 20,
                                           rmsGradientTolerance = 0.1 )


# ---------------------------------------------------------------------------
# (B) COMO O EasyHybrid CRIA O MODELO (trecho para session.py)
# ---------------------------------------------------------------------------
# Espelha exatamente o bloco do xTB/ORCA em define_a_new_QCModel. Basta um
# 'elif' novo:
#
#     elif parameters['qcengine'] == 'NWChem':
#         qcModel = QCModelNWChem.WithOptions (
#                       method     = parameters['method'    ],   # 'dft'
#                       functional = parameters['functional'],   # 'b3lyp'
#                       basis      = parameters['basis'      ],   # '6-31G*'
#                       keywords   = parameters['keywords'   ],   # linhas extra (lista) ou None
#                       memory     = parameters.get('memory', None),
#                       scratch    = parameters['scratch'    ],
#                   )
#         qcModel.randomScratch = parameters['random_scratch']
#         # arquivos a preservar por frame no backup (input/output):
#         system.e_nwchem_backup_files = parameters.get ( 'nwchem_backup_files', ['out'] )
#
# E no _common.py, registrar o backup no dispatcher (uma linha):
#     _QC_BACKUP_DISPATCH = {
#         'ORCA QC Model'   : _backup_orca_dispatch   ,
#         'XTB QC Model'    : _backup_xtb_dispatch    ,
#         'NWChem QC Model' : _backup_nwchem_dispatch ,   # <- novo
#     }


# ---------------------------------------------------------------------------
# (C) SISTEMA QC/MM  (regiao quantica + ambiente MM via point charges)
# ---------------------------------------------------------------------------
# Quando o sistema tem um nbModel (ambiente MM), o WriteInputFile emite
# automaticamente a linha 'bq load <arquivo> format 1 2 3 4 units bohr', e o
# NWChem inclui as cargas do ambiente no calculo QC. O fluxo e o mesmo:
#
#     system = ImportSystem ( "complexo.pkl" )         # sistema QC/MM ja montado
#     system.DefineQCModel ( qcModel, qcSelection = regiaoQC )
#     energy = system.Energy ( )
#
# [VALIDATE:QMMM] confirme que o arquivo de point charges que o pDynamo gera
# esta no layout 'x y z q' em Bohr. Se nao estiver, o WriteInputFile precisa
# converter antes de escrever a linha 'bq load'.
