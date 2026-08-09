#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: xtb_fragment_charges
#
#  Copyright 2022-2025 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
"""
xtb_fragment_charges
====================

Engine (Phase 1, NO graphical interface) to compute self-consistent per-fragment
MM charges using the `xtb` binary.

IDEA
----
The system is split into fragments (residue / segment / chain). In each CYCLE,
every fragment is computed in isolation with xTB, where:
  - covalent bonds cut at the fragment boundary receive a capping hydrogen
    (link atom), always appended AT THE END of the geometry;
  - the rest of the system enters as POINT CHARGES (electrostatic embedding),
    using the previous cycle's MM charges. In the 1st cycle, the other fragments
    are "mute" (no point charges).
The atomic charges (CM5 or Mulliken) are extracted, the boundary is treated, the
fragment is renormalized to its formal charge, and the result is written as MM
charges. This repeats until `max |Delta q| < tolerance`.

IMPORTANT SCIENTIFIC NOTE
-------------------------
CM5 is only available with GFN1 (GFN2 => Mulliken only). The engine validates
this: requesting CM5 with GFN2 emits a warning and falls back to Mulliken (or
forces GFN1, according to `cm5_policy`).

DEPENDS ON THE pDynamo SYSTEM only for: coordinates, connectivity
(Get12Indices) and initial MM charges (mmState.charges). Does NOT define a QCModel.

This module is testable in isolation (see `demo_dry_run` at the end) and supports
a DRY-RUN mode: it generates the xTB inputs and shows what it would do, without
executing the binary.
"""

import os
import re
import json
import math
import shutil
import tempfile
import subprocess


# --------------------------------------------------------------------------- #
#  Tabela minima de simbolos por numero atomico (fallback)                     #
# --------------------------------------------------------------------------- #

#_Z_TO_SYMBOL = {
#    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
#    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
#    16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 25: "Mn", 26: "Fe",
#    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 34: "Se", 35: "Br", 53: "I",
#}

_Z_TO_SYMBOL = {
     1: "H",   2: "He",
     3: "Li",  4: "Be",  5: "B",   6: "C",   7: "N",   8: "O",   9: "F",  10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",  16: "S",  17: "Cl", 18: "Ar",
    19: "K",  20: "Ca", 21: "Sc", 22: "Ti", 23: "V",  24: "Cr", 25: "Mn", 26: "Fe",
    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
    35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y",  40: "Zr", 41: "Nb", 42: "Mo",
    43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I",  54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce",
    59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy",
    67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf", 73: "Ta", 74: "W",
    75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb",
    83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
    91: "Pa", 92: "U",  93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf",
    99: "Es",100: "Fm",101: "Md",102: "No",103: "Lr",104: "Rf",105: "Db",106: "Sg",
   107: "Bh",108: "Hs",109: "Mt",110: "Ds",111: "Rg",112: "Cn",113: "Nh",114: "Fl",
   115: "Mc",116: "Lv",117: "Ts",118: "Og",
}



_RESIDUE_CHARGE = {
    "ALA":  0,
    "ARG": +1,
    "ASN":  0,
    "ASP": -1,
    "CYS":  0,
    "GLN":  0,
    "GLU": -1,
    "GLY":  0,
    "HIS":  0,   # see HID/HIE/HIP below
    "ILE":  0,
    "LEU":  0,
    "LYS": +1,
    "MET":  0,
    "PHE":  0,
    "PRO":  0,
    "SER":  0,
    "THR":  0,
    "TRP":  0,
    "TYR":  0,
    "VAL":  0,

    # Common protonation states
    "ASH":  0,   # protonated Asp
    "GLH":  0,   # protonated Glu
    "HID":  0,   # His delta-protonated
    "HIE":  0,   # His epsilon-protonated
    "HIP": +1,   # doubly protonated His
    "LYN":  0,   # neutral Lys
    "CYM": -1,   # deprotonated Cys
    "CYX":  0,   # disulfide Cys
    "TYM": -1,   # deprotonated Tyr

    "NTER": +1,   # NH3+
    "CTER": -1,   # COO-

}



ANGSTROM_TO_BOHR = 1.8897259886
# typical X-H distance for the capping hydrogen (Angstrom)
DEFAULT_CAP_DISTANCE = 1.09


# --------------------------------------------------------------------------- #
#  Estruturas de dados                                                         #
# --------------------------------------------------------------------------- #
class Fragment:
    """Um fragmento (residuo/segmento/cadeia) a ser calculado pelo xTB."""

    def __init__(self, key, atom_indexes):
        self.key = key                    # identificador (ex: 'A/ALA/12')
        self.atom_indexes = list(atom_indexes)  # indices GLOBAIS no system
        self.formal_charge = 0            # sugerido; editavel pelo usuario
        self.multiplicity = 1             # sugerido; editavel pelo usuario
        self.include = True               # if False, not recomputed (keeps MM charges)
        # filled in during the calculation:
        self.cap_atoms = []               # lista de dicts: H de capping
        # cap = {'pos': (x,y,z), 'caps_local_index': i, 'caps_global_index': g}
        self.last_charges = None          # atomic charges from the previous cycle

    def __repr__(self):
        return "Fragment({}, natoms={}, q={}, mult={})".format(
            self.key, len(self.atom_indexes), self.formal_charge, self.multiplicity)


# --------------------------------------------------------------------------- #
#  Acesso ao system pDynamo (isolado para facilitar teste/mock)               #
# --------------------------------------------------------------------------- #
class SystemAccessor:
    """Encapsula a leitura de dados do system pDynamo.

    Mantem toda a dependencia do pDynamo num so' lugar. Para testar sem
    pDynamo, basta fornecer um objeto com a mesma interface (ver o mock no
    bloco de demonstracao no fim do arquivo).
    """

    def __init__(self, system):
        self.system = system

    def natoms(self):
        return len(self.system.atoms)

    def symbol(self, index):
        atom = self.system.atoms[index]
        # pDynamo: atom.atomicNumber; fallbacks defensivos
        z = getattr(atom, "atomicNumber", None)
        if z is None:
            z = getattr(atom, "number", None)
        if z is not None:
            return _Z_TO_SYMBOL.get(int(z), "X")
        # ultimo recurso: primeira letra do label
        lab = getattr(atom, "label", "X")
        return lab[:1].upper() + (lab[1:2].lower() if len(lab) > 1 else "")

    def coordinates(self, index):
        """(x, y, z) em Angstrom para o atomo global `index`."""
        c = self.system.coordinates3[index]
        return (float(c[0]), float(c[1]), float(c[2]))

    def mm_charge(self, index):
        try:
            return float(self.system.mmState.charges[index])
        except Exception:
            return 0.0

    def set_mm_charge(self, index, value):
        self.system.mmState.charges[index] = float(value)

    def bonds_12(self):
        """Lista de pares (i, j) de atomos ligados (1-2), indices globais.

        Usa o mmState.mmTerms com Get12Indices() (mesma fonte usada no resto
        do EasyHybrid). Retorna [] se nao houver topologia MM.
        """
        pairs = []
        mm = getattr(self.system, "mmState", None)
        if mm is None:
            return pairs
        try:
            for term in mm.mmTerms:
                if hasattr(term, "Get12Indices"):
                    idx = list(term.Get12Indices())
                    # Get12Indices retorna uma lista achatada [i0,j0,i1,j1,...]
                    for k in range(0, len(idx) - 1, 2):
                        pairs.append((idx[k], idx[k + 1]))
                    break
        except Exception:
            pass
        return pairs


# --------------------------------------------------------------------------- #
#  Fragmentacao                                                                #
# --------------------------------------------------------------------------- #
def build_fragments(system, level="residue", residue_of=None, segment_of=None,
                    chain_of=None):
    """Divide o system em fragmentos e sugere carga formal e multiplicidade.

    Args:
        system: system pDynamo (ou mock com a mesma interface).
        level:  'residue' | 'segment' | 'chain'.
        residue_of/segment_of/chain_of: funcoes opcionais index->chave. Se None,
            tenta-se inferir do system (atom.residue etc.). Fornecer estas
            funcoes torna o modulo testavel sem pDynamo.

    Returns:
        lista de Fragment, com formal_charge sugerido = round(sum mm_charges) e
        multiplicity = 1.
    """
    #acc = SystemAccessor(system)
    #n = acc.natoms()
    #
    #key_fn = {"residue": residue_of, "segment": segment_of, "chain": chain_of}.get(level)
    #
    #if key_fn is None:
    #    key_fn = _default_key_fn(system, level)
    #
    #groups = {}
    #order = []
    #for i in range(n):
    #    k = key_fn(i)
    #    if k not in groups:
    #        groups[k] = []
    #        order.append(k)
    #    groups[k].append(i)
    #
    #fragments = []
    #for k in order:
    #    frag = Fragment(k, groups[k])
    #    print(frag)
    #    resn = frag[0]
    #    resn = resn.split('/')
    #    resn = resn[1]
    #    
    #    if resn in _RESIDUE_CHARGE.keys():
    #        qsum = _RESIDUE_CHARGE[resn]
    #        frag.formal_charge = int(round(qsum))
    #        
    #    else:
    #        # suggested formal charge = rounded sum of the MM charges
    #        qsum = sum(acc.mm_charge(i) for i in groups[k])
    #        frag.formal_charge = int(round(qsum))
    #    
    #    
    #    frag.multiplicity = 1
    #    fragments.append(frag)
    #return fragments
    acc = SystemAccessor(system)
    n = acc.natoms()

    key_fn = {
        "residue": residue_of,
        "segment": segment_of,
        "chain": chain_of
    }.get(level)

    if key_fn is None:
        key_fn = _default_key_fn(system, level)

    groups = {}
    order = []

    for i in range(n):
        key = key_fn(i)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(i)

    fragments = []

    for key in order:

        frag = Fragment(key, groups[key])

        # Residue name, if possible.
        residue_name = None
        if level == "residue":
            if isinstance(key, str):
                parts = key.split("/")
                residue_name = parts[-2]  # works for SEG/RES or just RES
            else:
                residue_name = str(key)
        #print(residue_name)
        
        if residue_name is not None and residue_name in _RESIDUE_CHARGE:
            frag.formal_charge = _RESIDUE_CHARGE[residue_name]
        else:
            qsum = sum(acc.mm_charge(i) for i in groups[key])
            frag.formal_charge = int(round(qsum))

        frag.multiplicity = 1
        fragments.append(frag)

    return fragments

def _default_key_fn(system, level):
    """Extrai a chave (residuo/segmento/cadeia) usando a API real do pDynamo3.

    No pDynamo3 a hierarquia NAO fica no atomo, e sim na sequencia:
      - atom.parent        -> o residuo/componente
      - atom.parent.parent -> a entity (cadeia/segmento); .label ex.: 'A' ou 'PRTA'
      - system.sequence.ParseLabel(atom.parent.label, fields=3) -> (resName, resSeq, iCode)

    Isto espelha _get_atom_info_from_pdynamo_atom_obj do EasyHybrid.
    """
    sequence = getattr(system, "sequence", None)

    def _chain_id(atom):
        try:
            return atom.parent.parent.label[0:1]
        except Exception:
            return "A"

    def _segment_id(atom):
        try:
            return atom.parent.parent.label[0:4]
        except Exception:
            return "SEG"

    def _residue_key(atom):
        # resName + resSeq via ParseLabel; cai para o label do parent se falhar
        try:
            if sequence is not None:
                resName, resSeq, iCode = sequence.ParseLabel(atom.parent.label, fields=3)
                return "{}/{}/{}".format(_chain_id(atom), resName, resSeq)
        except Exception:
            pass
        try:
            return "{}/{}".format(_chain_id(atom), atom.parent.label)
        except Exception:
            return "UNK"

    def fn(i):
        atom = system.atoms[i]
        if level == "chain":
            return _chain_id(atom)
        if level == "segment":
            return _segment_id(atom)
        return _residue_key(atom)

    return fn


# --------------------------------------------------------------------------- #
#  Capping (link atoms)                                                        #
# --------------------------------------------------------------------------- #
def add_caps(fragment, system, cap_distance=DEFAULT_CAP_DISTANCE):
    """Detecta ligacoes cortadas e adiciona H de capping AO FINAL do fragmento.

    Para cada ligacao (Q, M) onde Q pertence ao fragmento e M nao, coloca um H
    na direcao Q->M, a `cap_distance` Angstrom de Q. Guarda em fragment.cap_atoms
    de qual atomo pesado Q cada H faz o cap (necessario para o tratamento de
    carga na fronteira).
    """
    acc = SystemAccessor(system)
    in_frag = set(fragment.atom_indexes)
    bonds = acc.bonds_12()

    caps = []
    n_real = len(fragment.atom_indexes)
    for (a, b) in bonds:
        q, m = None, None
        if a in in_frag and b not in in_frag:
            q, m = a, b
        elif b in in_frag and a not in in_frag:
            q, m = b, a
        else:
            continue
        pq = acc.coordinates(q)
        pm = acc.coordinates(m)
        d = math.dist(pq, pm) if hasattr(math, "dist") else _dist(pq, pm)
        if d <= 1e-6:
            continue
        scale = cap_distance / d
        ph = (pq[0] + (pm[0] - pq[0]) * scale,
              pq[1] + (pm[1] - pq[1]) * scale,
              pq[2] + (pm[2] - pq[2]) * scale)
        # local H index = after all real atoms + caps already added
        caps.append({
            "pos": ph,
            "caps_global_index": q,                 # heavy atom (global) receiving the cap
            "caps_local_index": fragment.atom_indexes.index(q),  # idx local do Q
            "cap_local_index": n_real + len(caps),  # idx local deste H
        })
    fragment.cap_atoms = caps
    return caps


def _dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


# --------------------------------------------------------------------------- #
#  Geracao de input xTB                                                        #
# --------------------------------------------------------------------------- #
def write_xyz(fragment, system, path):
    """Escreve a geometria do fragmento (atomos reais + H de capping) em .xyz.

    Ordem: atomos reais (na ordem de fragment.atom_indexes), depois os H de
    capping (na ordem de fragment.cap_atoms). Assim os caps ficam SEMPRE AO FINAL.
    """
    acc = SystemAccessor(system)
    lines = []
    for gidx in fragment.atom_indexes:
        s = acc.symbol(gidx)
        x, y, z = acc.coordinates(gidx)
        lines.append("{:<3s} {:>15.8f} {:>15.8f} {:>15.8f}".format(s, x, y, z))
    for cap in fragment.cap_atoms:
        x, y, z = cap["pos"]
        lines.append("{:<3s} {:>15.8f} {:>15.8f} {:>15.8f}".format("H", x, y, z))

    natoms = len(fragment.atom_indexes) + len(fragment.cap_atoms)
    with open(path, "w") as fh:
        fh.write("{}\n".format(natoms))
        fh.write("fragment {}\n".format(fragment.key))
        fh.write("\n".join(lines) + "\n")
    return path


def write_point_charges(fragment, system, all_fragments, charges_by_index, path):
    """Escreve o arquivo de cargas pontuais (embedding) para o xTB.

    Formato xtb (`pcharge`, ativado por --input com $embedding, ou arquivo
    'pcharge'): primeira linha = numero de cargas; demais: q x y z [gam]
    com coordenadas em BOHR e carga em e.

    Inclui TODOS os atomos que NAO pertencem a este fragmento, usando
    `charges_by_index` (as cargas MM do ciclo anterior). Retorna o numero de
    cargas escritas (0 => nao gera embedding, ex.: 1o ciclo com tudo mudo).
    """
    acc = SystemAccessor(system)
    in_frag = set(fragment.atom_indexes)

    rows = []
    for gidx, q in charges_by_index.items():
        if gidx in in_frag:
            continue
        if abs(q) < 1e-9:
            continue
        x, y, z = acc.coordinates(gidx)
        rows.append((q, x * ANGSTROM_TO_BOHR, y * ANGSTROM_TO_BOHR, z * ANGSTROM_TO_BOHR))

    if not rows:
        return 0

    with open(path, "w") as fh:
        fh.write("{}\n".format(len(rows)))
        for (q, x, y, z) in rows:
            # o 5o campo (gam) e' opcional; omitido usa o default do xtb
            fh.write("{:>14.8f} {:>16.8f} {:>16.8f} {:>16.8f}\n".format(q, x, y, z))
    return len(rows)


def build_xtb_command(xtb_path, xyz_file, charge, multiplicity, method="gfn1",
                      charge_model="cm5", pcharge_file=None, extra=None):
    """Monta a lista de argumentos para chamar o binario xtb.

    method: 'gfn0' | 'gfn1' | 'gfn2'
    charge_model: 'cm5' | 'mulliken' (afeta so' a LEITURA; ver validacao)
    """
    gfn = {"gfn0": "0", "gfn1": "1", "gfn2": "2"}.get(method.lower(), "1")
    uhf = max(int(multiplicity) - 1, 0)   # numero de eletrons desemparelhados
    cmd = [xtb_path, os.path.basename(xyz_file),
           "--gfn", gfn,
           "--chrg", str(int(charge)),
           "--uhf", str(uhf),
           "--sp"]
    if pcharge_file:
        # xtb reads point charges from a file named 'pcharge' in the directory,
        # or via --input; here we assume the 'pcharge' file in the cwd.
        pass
    if extra:
        cmd += list(extra)
    return cmd


def validate_method_and_charge_model(method, charge_model, cm5_policy="fallback"):
    """CM5 so' existe em GFN1. Resolve o conflito conforme a politica.

    cm5_policy:
        'fallback'  -> se CM5 pedido com GFN2/0, usa Mulliken (com aviso).
        'force_gfn1'-> se CM5 pedido, forca method='gfn1'.
    Retorna (method, charge_model, aviso_ou_None).
    """
    cm = charge_model.lower()
    me = method.lower()
    if cm == "cm5" and me != "gfn1":
        if cm5_policy == "force_gfn1":
            return "gfn1", "cm5", "CM5 requer GFN1: metodo forcado para GFN1."
        return me, "mulliken", "CM5 indisponivel em {}: usando Mulliken.".format(me.upper())
    return me, cm, None


# --------------------------------------------------------------------------- #
#  Execucao e leitura                                                          #
# --------------------------------------------------------------------------- #
def run_xtb(cmd, workdir, dry_run=False):
    """Executa o xtb em `workdir`. Em dry_run, so' retorna o comando montado."""
    if dry_run:
        return {"dry_run": True, "cmd": " ".join(cmd), "workdir": workdir}
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    return {"returncode": proc.returncode, "stdout": proc.stdout, "stderr": proc.stderr}


def read_charges(workdir, charge_model="cm5", stdout_text=None):
    """Le as cargas atomicas do resultado do xtb.

    Prefere o `xtbout.json` (mais robusto). Se nao existir, tenta o
    `charges` (Mulliken, um valor por linha) e, para CM5, o texto do stdout
    ('Mulliken/CM5' na tabela de propriedades do GFN1).
    """
    # 1) JSON
    jpath = os.path.join(workdir, "xtbout.json")
    if os.path.isfile(jpath):
        try:
            with open(jpath) as fh:
                data = json.load(fh)
            # chaves possiveis conforme versao do xtb
            for key in ("partial charges", "mulliken partial charges", "charges"):
                if key in data:
                    return [float(x) for x in data[key]]
        except Exception:
            pass

    # 2) 'charges' file (Mulliken, 1 per line)
    cpath = os.path.join(workdir, "charges")
    if charge_model.lower() != "cm5" and os.path.isfile(cpath):
        try:
            with open(cpath) as fh:
                return [float(line.split()[0]) for line in fh if line.strip()]
        except Exception:
            pass

    # 3) stdout: tabela "Mulliken/CM5" do GFN1
    if stdout_text:
        return _parse_charges_from_stdout(stdout_text, charge_model)

    raise RuntimeError("Nao consegui ler as cargas do xtb em {}".format(workdir))


def _parse_charges_from_stdout(text, charge_model):
    """Extrai Mulliken ou CM5 da tabela do GFN1 no stdout.

    Formato tipico (GFN1):
        #   Z          Mulliken/CM5     ...
        1   O   ...    0.67569  0.33312  ...
    Mulliken = 1a coluna numerica, CM5 = 2a.
    """
    col = 1 if charge_model.lower() == "cm5" else 0
    charges = []
    started = False
    for line in text.splitlines():
        if "Mulliken/CM5" in line:
            started = True
            continue
        if started:
            parts = line.split()
            # linha valida: idx symbol num num ...
            nums = [p for p in parts if _is_float(p)]
            if len(nums) >= 2:
                charges.append(float(nums[col]))
            elif charges:
                break
    if not charges:
        raise RuntimeError("Tabela Mulliken/CM5 nao encontrada no stdout.")
    return charges


def _is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# --------------------------------------------------------------------------- #
#  Boundary treatment and renormalization                                      #
# --------------------------------------------------------------------------- #
def handle_boundary(fragment, charges, mode="redistribute"):
    """Trata as cargas dos H de capping (que nao sao atomos reais).

    charges: lista com natoms_reais + n_caps valores (caps ao final).
    mode:
      'redistribute' -> soma a carga de cada H de cap ao atomo pesado que ele
                        capeia (Q). Depois remove os caps. (opcao 1)
      'discard'      -> descarta a carga dos caps; a renormalizacao posterior
                        redistribui a diferenca por todos os atomos. (opcao 2)
      'keep_on_boundary' -> mantem a carga do cap no atomo de fronteira Q
                        (equivalente a 'redistribute' aqui, mas SEM renormalizar
                        depois pelo caller se assim desejado). (opcao 3)

    Retorna uma lista de cargas SO' para os atomos reais do fragmento
    (comprimento = len(fragment.atom_indexes)).
    """
    n_real = len(fragment.atom_indexes)
    real = list(charges[:n_real])
    caps = list(charges[n_real:])

    if mode in ("redistribute", "keep_on_boundary"):
        for cap, qcap in zip(fragment.cap_atoms, caps):
            li = cap["caps_local_index"]
            real[li] += qcap
    elif mode == "discard":
        pass  # simply ignore the cap charges
    else:
        raise ValueError("modo de fronteira desconhecido: {}".format(mode))
    return real


def renormalize(charges, formal_charge):
    """Ajusta as cargas para somar exatamente `formal_charge`.

    Distribui a diferenca igualmente por todos os atomos (esquema simples e
    estavel). Retorna nova lista.
    """
    n = len(charges)
    if n == 0:
        return charges
    diff = formal_charge - sum(charges)
    delta = diff / n
    return [q + delta for q in charges]


# --------------------------------------------------------------------------- #
#  Worker paralelo (nivel de modulo -> picklable para multiprocessing.Pool)     #
# --------------------------------------------------------------------------- #
def _run_one_fragment_worker(job):
    """Executa o xTB para UM fragmento, num processo separado.

    'job' e' um dict com SO' dados serializaveis (nada de objetos pDynamo):
        xtb_path, method, charge_model, charge, multiplicity,
        xyz_text, pcharge_text (ou None), workdir, n_real, cap_local_indices,
        cap_caps_local_indices, formal_charge, boundary_mode

    Retorna: (frag_pos, real_charges) onde real_charges ja' passou por
    tratamento de fronteira + renormalizacao. Em erro, retorna (frag_pos, None,
    mensagem).
    """
    try:
        os.makedirs(job["workdir"], exist_ok=True)
        xyz_path = os.path.join(job["workdir"], "mol.xyz")
        with open(xyz_path, "w") as fh:
            fh.write(job["xyz_text"])
        if job.get("pcharge_text"):
            with open(os.path.join(job["workdir"], "pcharge"), "w") as fh:
                fh.write(job["pcharge_text"])

        cmd = [job["xtb_path"], "mol.xyz",
               "--gfn", job["gfn"],
               "--chrg", str(int(job["charge"])),
               "--uhf", str(max(int(job["multiplicity"]) - 1, 0)),
               "--sp"]
        proc = subprocess.run(cmd, cwd=job["workdir"], capture_output=True, text=True)

        raw = read_charges(job["workdir"], charge_model=job["charge_model"],
                           stdout_text=proc.stdout)

        # boundary treatment (replicated here to avoid depending on the Fragment object)
        n_real = job["n_real"]
        real = list(raw[:n_real])
        caps = list(raw[n_real:])
        if job["boundary_mode"] in ("redistribute", "keep_on_boundary"):
            for li, qcap in zip(job["cap_caps_local_indices"], caps):
                real[li] += qcap
        # renormaliza
        diff = job["formal_charge"] - sum(real)
        delta = diff / len(real) if real else 0.0
        real = [q + delta for q in real]

        return (job["frag_pos"], real, None)
    except Exception as e:
        return (job.get("frag_pos"), None, "{}: {}".format(type(e).__name__, e))

# --------------------------------------------------------------------------- #
#  Self-consistent loop                                                        #
# --------------------------------------------------------------------------- #
def run_self_consistent(system, fragments, xtb_path,
                        method="gfn1", charge_model="cm5",
                        boundary_mode="redistribute",
                        tolerance=0.01, max_cycles=25,
                        cm5_policy="fallback",
                        dry_run=False, workroot=None, verbose=True,
                        progress_cb=None):
    """Executes the self-consistent cycle of MM charges by fragment.

    Returns:
        dict containing: 'converged' (bool), 'cycles' (int), 'max_dq' (float),
        'charges_by_index' (dict mapping global index -> final charge),
        and 'history'.
    """
    method, charge_model, warn = validate_method_and_charge_model(
        method, charge_model, cm5_policy)
    if warn and verbose:
        print("[warning]", warn)

    acc = SystemAccessor(system)
    n = acc.natoms()

    # Current MM charges (cycle 0). In the 1st cycle, other fragments are
    # MUTE => 0.
    charges_by_index = {i: 0.0 for i in range(n)}

    # Prepare caps only once (the geometry does not change between cycles).
    for frag in fragments:
        add_caps(frag, system)

    if workroot is None:
        workroot = tempfile.mkdtemp(prefix="xtb_fragq_")
    os.makedirs(workroot, exist_ok=True)

    history = []
    converged = False
    max_dq = None

    for cycle in range(1, max_cycles + 1):
        new_charges = dict(charges_by_index)
        cycle_dq = 0.0

        for fi, frag in enumerate(fragments):
            wdir = os.path.join(workroot, "cycle{:02d}_frag{:03d}".format(cycle, fi))
            os.makedirs(wdir, exist_ok=True)

            xyz = write_xyz(frag, system, os.path.join(wdir, "mol.xyz"))
            npc = write_point_charges(frag, system, fragments,
                                      charges_by_index if cycle > 1 else {},
                                      os.path.join(wdir, "pcharge"))
            cmd = build_xtb_command(
                xtb_path, xyz, frag.formal_charge, frag.multiplicity,
                method=method, charge_model=charge_model,
                pcharge_file=(os.path.join(wdir, "pcharge") if npc else None))

            res = run_xtb(cmd, wdir, dry_run=dry_run)

            if dry_run:
                history.append({"cycle": cycle, "fragment": frag.key,
                                "cmd": res["cmd"], "n_pointcharges": npc,
                                "n_caps": len(frag.cap_atoms)})
                continue

            raw = read_charges(wdir, charge_model=charge_model,
                               stdout_text=res.get("stdout"))
            real = handle_boundary(frag, raw, mode=boundary_mode)
            real = renormalize(real, frag.formal_charge)

            # Check convergence by comparing with the previous cycle for
            # this fragment.
            if frag.last_charges is not None:
                dq = max(abs(a - b) for a, b in zip(real, frag.last_charges))
                cycle_dq = max(cycle_dq, dq)
            frag.last_charges = real

            # Store the new MM charges for the real atoms in this fragment.
            for li, gidx in enumerate(frag.atom_indexes):
                new_charges[gidx] = real[li]

            if progress_cb:
                progress_cb(cycle, fi, frag, real)

        charges_by_index = new_charges

        if dry_run:
            # In dry-run mode, convergence is not evaluated; only one
            # demonstrative cycle is executed.
            return {"dry_run": True, "cycles": 1, "history": history,
                    "workroot": workroot}

        max_dq = cycle_dq
        history.append({"cycle": cycle, "max_dq": max_dq})
        if verbose:
            print("cycle {}: max |dq| = {:.5f}".format(cycle, max_dq))

        if cycle > 1 and max_dq is not None and max_dq < tolerance:
            converged = True
            break

    return {"converged": converged, "cycles": cycle, "max_dq": max_dq,
            "charges_by_index": charges_by_index, "history": history,
            "workroot": workroot}

## --------------------------------------------------------------------------- #
##  Loop auto-consistente                                                       #
## --------------------------------------------------------------------------- #
#def run_self_consistent(system, fragments, xtb_path,
#                        method="gfn1", charge_model="cm5",
#                        boundary_mode="redistribute",
#                        tolerance=0.01, max_cycles=25,
#                        cm5_policy="fallback",
#                        dry_run=False, workroot=None, verbose=True,
#                        progress_cb=None):
#    """Run the self-consistent per-fragment MM charge cycle.
#
#    Returns:
#        dict com: 'converged' (bool), 'cycles' (int), 'max_dq' (float),
#        'charges_by_index' (dict global index -> carga final), 'history'.
#    """
#    method, charge_model, warn = validate_method_and_charge_model(
#        method, charge_model, cm5_policy)
#    if warn and verbose:
#        print("[aviso]", warn)
#
#    acc = SystemAccessor(system)
#    n = acc.natoms()
#
#    # current MM charges (cycle 0). 1st cycle: other fragments MUTE => 0.
#    charges_by_index = {i: 0.0 for i in range(n)}
#
#    # prepare caps once (geometry does not change between cycles)
#    for frag in fragments:
#        add_caps(frag, system)
#
#    if workroot is None:
#        workroot = tempfile.mkdtemp(prefix="xtb_fragq_")
#    os.makedirs(workroot, exist_ok=True)
#
#    history = []
#    converged = False
#    max_dq = None
#
#    for cycle in range(1, max_cycles + 1):
#        new_charges = dict(charges_by_index)
#        cycle_dq = 0.0
#
#        for fi, frag in enumerate(fragments):
#            wdir = os.path.join(workroot, "cycle{:02d}_frag{:03d}".format(cycle, fi))
#            os.makedirs(wdir, exist_ok=True)
#
#            xyz = write_xyz(frag, system, os.path.join(wdir, "mol.xyz"))
#            npc = write_point_charges(frag, system, fragments,
#                                      charges_by_index if cycle > 1 else {},
#                                      os.path.join(wdir, "pcharge"))
#            cmd = build_xtb_command(
#                xtb_path, xyz, frag.formal_charge, frag.multiplicity,
#                method=method, charge_model=charge_model,
#                pcharge_file=(os.path.join(wdir, "pcharge") if npc else None))
#
#            res = run_xtb(cmd, wdir, dry_run=dry_run)
#
#            if dry_run:
#                history.append({"cycle": cycle, "fragment": frag.key,
#                                "cmd": res["cmd"], "n_pointcharges": npc,
#                                "n_caps": len(frag.cap_atoms)})
#                continue
#
#            raw = read_charges(wdir, charge_model=charge_model,
#                               stdout_text=res.get("stdout"))
#            real = handle_boundary(frag, raw, mode=boundary_mode)
#            real = renormalize(real, frag.formal_charge)
#
#            # measure convergence by comparing with this fragment's previous cycle
#            if frag.last_charges is not None:
#                dq = max(abs(a - b) for a, b in zip(real, frag.last_charges))
#                cycle_dq = max(cycle_dq, dq)
#            frag.last_charges = real
#
#            # store the new MM charges of this fragment's real atoms
#            for li, gidx in enumerate(frag.atom_indexes):
#                new_charges[gidx] = real[li]
#
#            if progress_cb:
#                progress_cb(cycle, fi, frag, real)
#
#        charges_by_index = new_charges
#
#        if dry_run:
#            # in dry-run there is no convergence; runs only 1 demo cycle
#            return {"dry_run": True, "cycles": 1, "history": history,
#                    "workroot": workroot}
#
#        max_dq = cycle_dq
#        history.append({"cycle": cycle, "max_dq": max_dq})
#        if verbose:
#            print("cycle {}: max |dq| = {:.5f}".format(cycle, max_dq))
#
#        if cycle > 1 and max_dq is not None and max_dq < tolerance:
#            converged = True
#            break
#
#    return {"converged": converged, "cycles": cycle, "max_dq": max_dq,
#            "charges_by_index": charges_by_index, "history": history,
#            "workroot": workroot}



def run_self_consistent_parallel(system, fragments, xtb_path,
                                 method="gfn1", charge_model="cm5",
                                 boundary_mode="redistribute",
                                 tolerance=0.01, max_cycles=25,
                                 cm5_policy="fallback", nprocs=4,
                                 workroot=None, verbose=True,
                                 cycle_cb=None, init_from_mm=True):
    """PARALLEL version of the loop: within each cycle, fragments are run in
    separate processes (multiprocessing.Pool), since they are independent
    (all read the charges from the previous cycle). There is a barrier
    between cycles.

    Fragments with frag.include == False are NOT recalculated: they retain
    their current charges (by default, the original MM charges from the
    system) and still CONTRIBUTE to the electrostatic embedding of the other
    fragments. This prevents creating an electrostatic "hole" where the
    ignored fragment is located.

    init_from_mm: if True, the embedding in cycle 1 uses the ORIGINAL MM
    charges of the ignored fragments (recommended). If False, everything
    starts at zero.

    cycle_cb(cycle, max_dq): optional callback called at the end of each cycle.
    Returns the same dict as run_self_consistent.
    """
    import multiprocessing

    method, charge_model, warn = validate_method_and_charge_model(
        method, charge_model, cm5_policy)
    if warn and verbose:
        print("[warning]", warn)
    gfn = {"gfn0": "0", "gfn1": "1", "gfn2": "2"}.get(method.lower(), "1")

    acc = SystemAccessor(system)
    n = acc.natoms()

    # Current charges. Ignored fragments retain their original MM charge;
    # fragments that will be calculated start at 0 (neutral in cycle 1).
    included_indexes = set()
    for frag in fragments:
        if getattr(frag, "include", True):
            included_indexes.update(frag.atom_indexes)

    charges_by_index = {}
    for i in range(n):
        if i in included_indexes:
            charges_by_index[i] = 0.0
        else:
            # Atom belonging to an ignored fragment (or outside any fragment):
            # retain the original MM charge to contribute to the embedding.
            charges_by_index[i] = acc.mm_charge(i) if init_from_mm else 0.0

    # Only process fragments marked for inclusion.
    active_fragments = [f for f in fragments if getattr(f, "include", True)]

    for frag in active_fragments:
        add_caps(frag, system)

    if workroot is None:
        workroot = tempfile.mkdtemp(prefix="xtb_fragq_")
    os.makedirs(workroot, exist_ok=True)

    # Pre-extract the geometry of each fragment ONCE (it does not change
    # between cycles).
    frag_xyz_text = []
    frag_meta = []
    for frag in active_fragments:
        lines = []
        for gidx in frag.atom_indexes:
            s = acc.symbol(gidx)
            x, y, z = acc.coordinates(gidx)
            lines.append("{:<3s} {:>15.8f} {:>15.8f} {:>15.8f}".format(s, x, y, z))
        for cap in frag.cap_atoms:
            x, y, z = cap["pos"]
            lines.append("{:<3s} {:>15.8f} {:>15.8f} {:>15.8f}".format("H", x, y, z))
        natoms = len(frag.atom_indexes) + len(frag.cap_atoms)
        xyz_text = "{}\nfragment {}\n{}\n".format(natoms, frag.key, "\n".join(lines))
        frag_xyz_text.append(xyz_text)
        frag_meta.append({
            "n_real": len(frag.atom_indexes),
            "cap_caps_local_indices": [c["caps_local_index"] for c in frag.cap_atoms],
        })

    history = []
    converged = False
    max_dq = None
    last_frag_charges = [None] * len(active_fragments)

    for cycle in range(1, max_cycles + 1):
        jobs = []
        for fpos, frag in enumerate(active_fragments):
            wdir = os.path.join(workroot, "cycle{:02d}_frag{:03d}".format(cycle, fpos))
            pcharge_text = None
            if cycle > 1:
                in_frag = set(frag.atom_indexes)
                rows = []
                for gidx, q in charges_by_index.items():
                    if gidx in in_frag or abs(q) < 1e-9:
                        continue
                    x, y, z = acc.coordinates(gidx)
                    rows.append("{:>14.8f} {:>16.8f} {:>16.8f} {:>16.8f}".format(
                        q, x * ANGSTROM_TO_BOHR, y * ANGSTROM_TO_BOHR, z * ANGSTROM_TO_BOHR))
                if rows:
                    pcharge_text = "{}\n{}\n".format(len(rows), "\n".join(rows))

            jobs.append({
                "frag_pos": fpos, "workdir": wdir, "xtb_path": xtb_path,
                "gfn": gfn, "charge_model": charge_model,
                "charge": frag.formal_charge, "multiplicity": frag.multiplicity,
                "xyz_text": frag_xyz_text[fpos], "pcharge_text": pcharge_text,
                "n_real": frag_meta[fpos]["n_real"],
                "cap_caps_local_indices": frag_meta[fpos]["cap_caps_local_indices"],
                "formal_charge": frag.formal_charge,
                "boundary_mode": boundary_mode,
            })

        pool = multiprocessing.Pool(processes=nprocs)
        try:
            results = pool.map(_run_one_fragment_worker, jobs)
        finally:
            pool.close()
            pool.join()

        cycle_dq = 0.0
        new_charges = dict(charges_by_index)
        for (fpos, real, err) in results:
            if err:
                raise RuntimeError("fragment {} failed: {}".format(
                    active_fragments[fpos].key, err))
            if last_frag_charges[fpos] is not None:
                dq = max(abs(a - b) for a, b in zip(real, last_frag_charges[fpos]))
                cycle_dq = max(cycle_dq, dq)
            last_frag_charges[fpos] = real
            for li, gidx in enumerate(active_fragments[fpos].atom_indexes):
                new_charges[gidx] = real[li]

        charges_by_index = new_charges
        max_dq = cycle_dq
        history.append({"cycle": cycle, "max_dq": max_dq})
        if verbose:
            print("cycle {}: max |dq| = {:.5f}".format(cycle, max_dq))
        if cycle_cb:
            cycle_cb(cycle, max_dq)

        if cycle > 1 and max_dq < tolerance:
            converged = True
            break

    return {"converged": converged, "cycles": cycle, "max_dq": max_dq,
            "charges_by_index": charges_by_index, "history": history,
            "workroot": workroot}



#def run_self_consistent_parallel(system, fragments, xtb_path,
#                                 method="gfn1", charge_model="cm5",
#                                 boundary_mode="redistribute",
#                                 tolerance=0.01, max_cycles=25,
#                                 cm5_policy="fallback", nprocs=4,
#                                 workroot=None, verbose=True,
#                                 cycle_cb=None, init_from_mm=True):
#    """PARALLEL version of the loop: within each cycle, fragments run in
#    processos separados (multiprocessing.Pool), pois sao independentes (todos
#    read the previous cycle's charges). There is a barrier between cycles.
#
#    Fragmentos com frag.include == False NAO sao recalculados: eles mantem as
#    current charges (by default the system's original MM charges) and still
#    CONTRIBUTE to the electrostatic embedding of the others. This avoids creating a
#    electrostatic "hole" where the ignored fragment is.
#
#    init_from_mm: if True, cycle-1 embedding uses the ORIGINAL MM charges of the
#    fragmentos ignorados (recomendado). Se False, tudo comeca em zero.
#
#    cycle_cb(cycle, max_dq): optional callback at the end of each cycle.
#    Retorna o mesmo dict de run_self_consistent.
#    """
#    import multiprocessing
#
#    method, charge_model, warn = validate_method_and_charge_model(
#        method, charge_model, cm5_policy)
#    if warn and verbose:
#        print("[aviso]", warn)
#    gfn = {"gfn0": "0", "gfn1": "1", "gfn2": "2"}.get(method.lower(), "1")
#
#    acc = SystemAccessor(system)
#    n = acc.natoms()
#
#    # Cargas correntes. Fragmentos ignorados mantem sua carga MM original; os
#    # to be computed start at 0 (mute in cycle 1).
#    included_indexes = set()
#    for frag in fragments:
#        if getattr(frag, "include", True):
#            included_indexes.update(frag.atom_indexes)
#
#    charges_by_index = {}
#    for i in range(n):
#        if i in included_indexes:
#            charges_by_index[i] = 0.0
#        else:
#            # atom of an ignored fragment (or outside any fragment):
#            # keeps the original MM charge to take part in the embedding.
#            charges_by_index[i] = acc.mm_charge(i) if init_from_mm else 0.0
#
#    # so' processa os fragmentos marcados
#    active_fragments = [f for f in fragments if getattr(f, "include", True)]
#
#    for frag in active_fragments:
#        add_caps(frag, system)
#
#    if workroot is None:
#        workroot = tempfile.mkdtemp(prefix="xtb_fragq_")
#    os.makedirs(workroot, exist_ok=True)
#
#    # pre-extract each fragment's geometry ONCE (does not change between cycles)
#    frag_xyz_text = []
#    frag_meta = []
#    for frag in active_fragments:
#        lines = []
#        for gidx in frag.atom_indexes:
#            s = acc.symbol(gidx)
#            x, y, z = acc.coordinates(gidx)
#            lines.append("{:<3s} {:>15.8f} {:>15.8f} {:>15.8f}".format(s, x, y, z))
#        for cap in frag.cap_atoms:
#            x, y, z = cap["pos"]
#            lines.append("{:<3s} {:>15.8f} {:>15.8f} {:>15.8f}".format("H", x, y, z))
#        natoms = len(frag.atom_indexes) + len(frag.cap_atoms)
#        xyz_text = "{}\nfragment {}\n{}\n".format(natoms, frag.key, "\n".join(lines))
#        frag_xyz_text.append(xyz_text)
#        frag_meta.append({
#            "n_real": len(frag.atom_indexes),
#            "cap_caps_local_indices": [c["caps_local_index"] for c in frag.cap_atoms],
#        })
#
#    history = []
#    converged = False
#    max_dq = None
#    last_frag_charges = [None] * len(active_fragments)
#
#    for cycle in range(1, max_cycles + 1):
#        jobs = []
#        for fpos, frag in enumerate(active_fragments):
#            wdir = os.path.join(workroot, "cycle{:02d}_frag{:03d}".format(cycle, fpos))
#            pcharge_text = None
#            if cycle > 1:
#                in_frag = set(frag.atom_indexes)
#                rows = []
#                for gidx, q in charges_by_index.items():
#                    if gidx in in_frag or abs(q) < 1e-9:
#                        continue
#                    x, y, z = acc.coordinates(gidx)
#                    rows.append("{:>14.8f} {:>16.8f} {:>16.8f} {:>16.8f}".format(
#                        q, x * ANGSTROM_TO_BOHR, y * ANGSTROM_TO_BOHR, z * ANGSTROM_TO_BOHR))
#                if rows:
#                    pcharge_text = "{}\n{}\n".format(len(rows), "\n".join(rows))
#
#            jobs.append({
#                "frag_pos": fpos, "workdir": wdir, "xtb_path": xtb_path,
#                "gfn": gfn, "charge_model": charge_model,
#                "charge": frag.formal_charge, "multiplicity": frag.multiplicity,
#                "xyz_text": frag_xyz_text[fpos], "pcharge_text": pcharge_text,
#                "n_real": frag_meta[fpos]["n_real"],
#                "cap_caps_local_indices": frag_meta[fpos]["cap_caps_local_indices"],
#                "formal_charge": frag.formal_charge,
#                "boundary_mode": boundary_mode,
#            })
#
#        pool = multiprocessing.Pool(processes=nprocs)
#        try:
#            results = pool.map(_run_one_fragment_worker, jobs)
#        finally:
#            pool.close()
#            pool.join()
#
#        cycle_dq = 0.0
#        new_charges = dict(charges_by_index)
#        for (fpos, real, err) in results:
#            if err:
#                raise RuntimeError("fragment {} failed: {}".format(
#                    active_fragments[fpos].key, err))
#            if last_frag_charges[fpos] is not None:
#                dq = max(abs(a - b) for a, b in zip(real, last_frag_charges[fpos]))
#                cycle_dq = max(cycle_dq, dq)
#            last_frag_charges[fpos] = real
#            for li, gidx in enumerate(active_fragments[fpos].atom_indexes):
#                new_charges[gidx] = real[li]
#
#        charges_by_index = new_charges
#        max_dq = cycle_dq
#        history.append({"cycle": cycle, "max_dq": max_dq})
#        if verbose:
#            print("cycle {}: max |dq| = {:.5f}".format(cycle, max_dq))
#        if cycle_cb:
#            cycle_cb(cycle, max_dq)
#
#        if cycle > 1 and max_dq < tolerance:
#            converged = True
#            break
#
#    return {"converged": converged, "cycles": cycle, "max_dq": max_dq,
#            "charges_by_index": charges_by_index, "history": history,
#            "workroot": workroot}
#

def apply_charges_to_system(system = None, charges_by_index = None, factor = 1.0):
    """Writes the final charges back to the pDynamo system (mmState.charges)."""
    acc = SystemAccessor(system)
    for gidx, q in charges_by_index.items():
        acc.set_mm_charge(gidx, q*factor)


# --------------------------------------------------------------------------- #
#  Isolated demonstration/test (mock, without pDynamo or xTB)
# --------------------------------------------------------------------------- #
def _demo_mock_system():
    """Creates a minimal mock system that mimics the real pDynamo3 API.

    Structure: atom.parent = residue (with .label 'HOH.1'), 
    atom.parent.parent = entity/chain (with .label 'A'), 
    and system.sequence.ParseLabel(...) returns (resName, resSeq, iCode). 
    Two “water” molecules with a fictitious bond between them are 
    included to generate a capping boundary..
    """
    class _Entity:
        def __init__(self, label):
            self.label = label

    class _Residue:
        def __init__(self, label, entity):
            self.label = label
            self.parent = entity

    class _Atom:
        def __init__(self, z, residue):
            self.atomicNumber = z
            self.parent = residue
            self.label = "X"

    class _Sequence:
        def ParseLabel(self, label, fields=3):
            # 'HOH.1' -> ('HOH', '1', '')
            parts = label.split(".")
            resName = parts[0]
            resSeq = parts[1] if len(parts) > 1 else "0"
            return resName, resSeq, ""

    class _MMState:
        def __init__(self, charges, terms):
            self.charges = charges
            self.mmTerms = terms

    class _Term:
        def __init__(self, pairs):
            self._pairs = pairs
        def Get12Indices(self):
            flat = []
            for (i, j) in self._pairs:
                flat += [i, j]
            return flat

    class _Sys:
        pass

    chainA = _Entity("A")
    res1 = _Residue("HOH.1", chainA)
    res2 = _Residue("HOH.2", chainA)

    s = _Sys()
    s.atoms = [
        _Atom(8, res1), _Atom(1, res1), _Atom(1, res1),
        _Atom(8, res2), _Atom(1, res2), _Atom(1, res2),
    ]
    s.coordinates3 = [
        (0.000, 0.000, 0.000), (0.757, 0.586, 0.000), (-0.757, 0.586, 0.000),
        (3.000, 0.000, 0.000), (3.757, 0.586, 0.000), (2.243, 0.586, 0.000),
    ]
    s.sequence = _Sequence()
    charges = [-0.8, 0.4, 0.4, -0.8, 0.4, 0.4]
    terms = [_Term([(0, 1), (0, 2), (3, 4), (3, 5), (2, 3)])]
    s.mmState = _MMState(charges, terms)
    return s


def demo_dry_run():
    """Roda o pipeline em DRY-RUN sobre o mock e imprime o que faria."""
    system = _demo_mock_system()
    frags = build_fragments(system, level="residue")
    print("Fragmentos sugeridos:")
    for f in frags:
        print("  ", f)
    result = run_self_consistent(
        system, frags, xtb_path="xtb",
        method="gfn1", charge_model="cm5",
        boundary_mode="redistribute", tolerance=0.01,
        dry_run=True, verbose=True)
    print("\nDRY-RUN — commands that would be executed:")
    for h in result["history"]:
        print("  [{}] {}  (caps={}, pointcharges={})".format(
            h["fragment"], h["cmd"], h["n_caps"], h["n_pointcharges"]))
    print("\nworkroot:", result["workroot"])


if __name__ == "__main__":
    demo_dry_run()
