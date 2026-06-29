#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Copyright 2022-2025 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      Provides functions for selecting atoms and residues in pDynamo systems
#      to facilitate QM/MM partitioning and molecular simulations.
#
from pprint import pprint
import re


def get_file_type(filename):
    file_type = filename.split('.')
    return file_type[-1]
    
def read_SIMPLE_txt (filein = None):
    """ Function doc """    
    charges = []
    
    with open(filein, 'r') as _file:
        filetext = _file.readlines()
        
        for line in filetext:
            #print(line)
            #line = line.strip()
            if len(line) > 0:
                try:
                    line = float(line)
                    charges.append(line)
                except:
                    charges = False
                    #break
            else:
                pass
    return charges
    
def read_MOPAC_out (filein = None):
    """ Function doc """

def read_MOPAC_aux (filein = None):
    """ Function doc """
    data = open(filein, 'r')
    data = data.read()
    
    #.ATOM_CORE
    try:
        data2 = data.split('ATOM_CORE')
        data2 = data2[-1].split('ATOM_X:ANGSTROMS')
        data2 = data2[0] 

        data2 = data2.split('=')
        data2 = data2[-1].split()

        ATOM_CORE = data2
    except:
        ATOM_CORE = None

    #.HEAT_OF_FORMATION
    try:
        data2 = data.split('HEAT_OF_FORMATION:KCAL/MOL')
        data2 = data2[-1].split('GRADIENT_NORM:KCAL/MOL/ANGSTROM')
        data2 = data2[0].replace('=','')
        data2 = data2.replace('D','E')
        data2 = data2.replace('+','')
        data2 = data2.replace('\n','')
        data2 = float(data2)

        HEAT_OF_FORMATION = data2
    except:
        HEAT_OF_FORMATION = None


    #.ATOM_X_OPT
    try:
        data2 = data.split('ATOM_X_OPT:ANGSTROMS')
        data2 = data2[-1].split('ATOM_CHARGES')
        data2 = data2[0] 

        data2 = data2.split('=')
        data2 = data2[-1].split()

        ATOM_X_OPT = data2
    except:
        ATOM_X_OPT = None
        
        
    #,GRADIENTS
    try: 
        data2 = data.split('GRADIENTS:KCAL/MOL/ANGSTROM')
        data2 = data2[-1].split('OVERLAP_MATRIX')
        data2 = data2[0] 
        data2 = data2.split('=')
        data2 = data2[-1].split()

        GRADIENTS = data2
    except:
        GRADIENTS = None
    
    #,charges
    try: 
        data2 = data.split('ATOM_CHARGES')
        data2 = data2[-1].split('AO_CHARGES')
        data2 = data2[0] 
        data2 = data2.split('=')
        data2 = data2[-1].split()
                
        CHARGES = data2
    except:
        CHARGES = None    
    
    #print (HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS)
    return HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS, CHARGES
    
def read_MOPAC_arc (filein = None):
    """ Function doc """
    
def read_ORCA_log (logfile = None):
    """
    Extrai cargas atômicas de arquivos ORCA:
      - Mulliken
      - Loewdin
      - CHELPG
      - RESP (se disponível)

    Returns
    -------
    dict
        {
            "mulliken": [(idx, atom, charge), ...],
            "loewdin": [...],
            "chelpg": [...],
            "resp": [...]
        }
    """

    patterns = {
        "mulliken": re.compile(r"MULLIKEN ATOMIC CHARGES", re.I),
        "loewdin": re.compile(r"LOEWDIN ATOMIC CHARGES", re.I),
        "chelpg": re.compile(r"CHELPG Charges", re.I),
        "resp": re.compile(r"RESP Charges", re.I),
    }

    charges = {k: [] for k in patterns}

    with open(logfile, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]

        matched_key = None
        for key, pat in patterns.items():
            if pat.search(line):
                matched_key = key
                break

        if matched_key:
            i += 1

            # pula cabeçalhos / linhas vazias
            while i < len(lines):
                if lines[i].strip() == "":
                    i += 1
                    continue
                if set(lines[i].strip()) <= {"-", "="}:
                    i += 1
                    continue
                break

            # lê bloco
            while i < len(lines):
                row = lines[i].strip()

                if not row:
                    break

                # Caso Mulliken/Loewdin:
                # 0 C :  -0.12345
                m = re.match(
                    r"^\s*(\d+)\s+([A-Za-z]+)\s*[: ]+\s*(-?\d+\.\d+)",
                    row
                )

                # Caso CHELPG:
                # 0 C -0.123
                if not m:
                    m = re.match(
                        r"^\s*(\d+)\s+([A-Za-z]+)\s+(-?\d+\.\d+)",
                        row
                    )

                if not m:
                    break

                idx = int(m.group(1))
                atom = m.group(2)
                charge = float(m.group(3))

                #charges[matched_key].append((idx, atom, charge))
                charges[matched_key].append( charge )
                i += 1

        i += 1

    return charges

def read_XTB_log (logfile = None):
#def extract_xtb_charges(logfile):
    """
    Extrai cargas do output do xTB.

    Retorna:
    {
        "mulliken": [(idx, Z, charge), ...],
        "cm5": [(idx, atom, charge), ...]
    }
    """

    data = {
        "mulliken": [],
        "cm5": []
    }

    with open(logfile, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]

        # -------------------------
        # Mulliken block
        # -------------------------
        if "Mulliken population" in line:
            i += 1

            # pula header
            while i < n:
                if lines[i].strip().startswith("#"):
                    i += 1
                    break
                i += 1

            while i < n:
                row = lines[i].strip()

                if not row:
                    break

                parts = row.split()

                # Esperado:
                # idx Z covCN q ...
                if len(parts) >= 4:
                    try:
                        idx = int(parts[0])
                        atomic_number = int(parts[1])
                        charge = float(parts[3])

                        data["mulliken"].append(
                            #(idx, atomic_number, charge)
                             charge
                        )
                    except ValueError:
                        break
                else:
                    break

                i += 1

        # -------------------------
        # CM5 block
        # -------------------------
        elif "CM5 charges" in line:
            i += 1

            while i < n:
                row = lines[i].strip()

                if not row:
                    break

                m = re.match(
                    r"^\s*(\d+)\s+([A-Za-z]+)\s+(-?\d+\.\d+)",
                    row
                )

                if not m:
                    break

                idx = int(m.group(1))
                atom = m.group(2)
                charge = float(m.group(3))

                #data["cm5"].append((idx, atom, charge))
                data["cm5"].append(charge)
                i += 1

        i += 1

    return data

def read_MOL2 (filein):
    """ Function doc """
    #at  = MolecularProperties.AtomTypes()
    
    with open(filein, 'r') as mol2_file:
        filetext = mol2_file.read()

        molecules     =  filetext.split('@<TRIPOS>MOLECULE')
        firstmolecule =  molecules[1].split('@<TRIPOS>ATOM')
        header        =  firstmolecule[0]
        firstmolecule =  firstmolecule[1].split('@<TRIPOS>BOND')
        raw_atoms     =  firstmolecule[0]
        bonds         =  firstmolecule[1]


    header    = header.split('\n')
    raw_atoms = raw_atoms.split('\n')
    bonds     = bonds.split('\n')
    
    atoms = {
            'id'     : [],
            'name'   : [],
            'xyz'    : [],
            'type'   : [],
            'resi'   : [],
            'resn'   : [],
            'charges': [],
            }
    
    _bonds = []
    
    #pprint (raw_atoms)
    for atom in raw_atoms[1:]:
        if len(atom) == 0:
            pass
        else:
            try:
                a = atom.split()
                atoms['id'     ].append(int(a[0]))
                atoms['name'   ].append(a[1])
                atoms['xyz'    ].append([ float(a[2]), float(a[3]), float(a[4])] )
                atoms['type'   ].append(a[5])
                atoms['resi'   ].append(int(a[6]))
                atoms['resn'   ].append(a[7])
                atoms['charges'].append(float(a[8]))
            except:
                pass
    #pprint (atoms)
    
    
    for bond in bonds:
        b = bond.split()
        
        if len(b) == 4:
            _bonds.append(b)
            
        else:
            pass
    #pprint (_bonds)
    return atoms, _bonds


def extract_orca_charges (filename):
    """
    Extrai cargas atômicas de arquivos ORCA:
      - Mulliken
      - Loewdin
      - CHELPG
      - RESP (se disponível)

    Returns
    -------
    dict
        {
            "mulliken": [(idx, atom, charge), ...],
            "loewdin": [...],
            "chelpg": [...],
            "resp": [...]
        }
    """

    patterns = {
        "mulliken": re.compile(r"MULLIKEN ATOMIC CHARGES", re.I),
        "loewdin": re.compile(r"LOEWDIN ATOMIC CHARGES", re.I),
        "chelpg": re.compile(r"CHELPG Charges", re.I),
        "resp": re.compile(r"RESP Charges", re.I),
    }

    charges = {k: [] for k in patterns}

    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]

        matched_key = None
        for key, pat in patterns.items():
            if pat.search(line):
                matched_key = key
                break

        if matched_key:
            i += 1

            # pula cabeçalhos / linhas vazias
            while i < len(lines):
                if lines[i].strip() == "":
                    i += 1
                    continue
                if set(lines[i].strip()) <= {"-", "="}:
                    i += 1
                    continue
                break

            # lê bloco
            while i < len(lines):
                row = lines[i].strip()

                if not row:
                    break

                # Caso Mulliken/Loewdin:
                # 0 C :  -0.12345
                m = re.match(
                    r"^\s*(\d+)\s+([A-Za-z]+)\s*[: ]+\s*(-?\d+\.\d+)",
                    row
                )

                # Caso CHELPG:
                # 0 C -0.123
                if not m:
                    m = re.match(
                        r"^\s*(\d+)\s+([A-Za-z]+)\s+(-?\d+\.\d+)",
                        row
                    )

                if not m:
                    break

                idx = int(m.group(1))
                atom = m.group(2)
                charge = float(m.group(3))

                #charges[matched_key].append((idx, atom, charge))
                charges[matched_key].append( charge )
                i += 1

        i += 1

    return charges

def extract_txt_charges (filename):
    """ Function doc """
    charges = []
    
    with open(filename, 'r') as _file:
        filetext = _file.readlines()
        
        for line in filetext:
            print(line)
            #line = line.strip()
            if len(line) > 0:
                try:
                    line = float(line)
                    charges.append(line)
                except:
                    charges = False
                    #break
            else:
                pass
    charges = {'unk': charges}
    return charges

def extract_xtb_charges(filename):
    """
    Extrai cargas Mulliken e CM5 de outputs recentes do xTB.

    Returns
    -------
    dict
        {
            "mulliken": [(idx, atom, charge)],
            "cm5": [(idx, atom, charge)],
            "orbital_populations": [
                (idx, atom, ns, np, nd)
            ]
        }
    """

    results = {
        "mulliken": [],
        "cm5": [],
        #"orbital_populations": []
    }

    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    inside_block = False

    pattern = re.compile(
        r"^\s*(\d+)([A-Za-z]+)\s+"      # 1C / 16O / 25Cl
        r"(-?\d+\.\d+)\s+"              # Mulliken
        r"(-?\d+\.\d+)\s+"              # CM5
        r"(\d+\.\d+)\s+"                # n(s)
        r"(\d+\.\d+)\s+"                # n(p)
        r"(\d+\.\d+)"                   # n(d)
    )

    for line in lines:

        if "Mulliken/CM5 charges" in line:
            inside_block = True
            continue

        if inside_block:
            stripped = line.strip()

            if not stripped:
                break

            match = pattern.match(line)

            if not match:
                break

            idx = int(match.group(1))
            atom = match.group(2)

            mulliken = float(match.group(3))
            cm5 = float(match.group(4))

            ns = float(match.group(5))
            np = float(match.group(6))
            nd = float(match.group(7))

            #results["mulliken"].append(
            #    (idx, atom, mulliken)
            results["mulliken"].append(
                mulliken
            )

            #results["cm5"].append(
            #    (idx, atom, cm5)
            results["cm5"].append(
                cm5
            )

            #results["orbital_populations"].append(
            #    (idx, atom, ns, np, nd)
            #)

    return results

def extract_aux_charges(filename):
    HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS, CHARGES = read_MOPAC_aux (filein = filename)
    
    charges = {'mulliken': CHARGES}
    return(charges)
    
def extract_mol2_charges(filename):
    atoms, _bonds = read_MOL2 (filename)
    charges = []
    for atom in atoms:
        charges.append(atoms['charges'])
    
    return charges
    
def extract_mopac_charges(filename):
    """
    Parser genérico para MOPAC (inclui versões 2016–23).

    Returns
    -------
    list[dict]
        [
            {
                "index": 1,
                "symbol": "C",
                "charge": -0.1234
            }
        ]
    """

    headers = [
        "NET ATOMIC CHARGES",
        "NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS",
        "ATOM NO.   TYPE          CHARGE",
        "ATOM   CHEMICAL SYMBOL   NET CHARGE",
    ]

    patterns = [
        # 1   C   -0.1234
        re.compile(
            r"^\s*(\d+)\s+([A-Za-z]{1,3})\s+(-?\d+\.\d+)"
        ),

        # atom=1 element=C charge=-0.1234
        re.compile(
            r".*?(\d+).*?([A-Za-z]{1,3}).*?(-?\d+\.\d+)"
        )
    ]

    atoms = []
    inside_block = False
    charges = []
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            stripped = line.strip()

            if any(header in line for header in headers):
                inside_block = True
                continue

            if not inside_block:
                continue

            if not stripped:
                if atoms:
                    break
                continue

            if set(stripped) <= {"-", "="}:
                continue

            matched = False

            for pattern in patterns:
                m = pattern.match(line)

                if m:
                    idx = int(m.group(1))
                    symbol = m.group(2)
                    charge = float(m.group(3))

                    atoms.append({
                        "index": idx,
                        "symbol": symbol,
                        "charge": charge
                    })
                    charges.append(charge)
                    matched = True
                    break

            if not matched and atoms:
                break
    
    charges = {'mulliken': charges}
    #return atoms, charges
    return charges

def detect_qm_log_type(filename):
    """
    Detecta se um arquivo de log pertence ao:
    - ORCA
    - xTB
    - MOPAC

    Returns
    -------
    str
        'orca', 'xtb', 'mopac', or 'unknown'
    """

    signatures = {
        "xtb": [
            "xTB version",
            "GFN2-xTB",
            "Mulliken/CM5 charges",
            "Wiberg/Mayer (AO) data",
            "covCN"
        ],

        "mopac": [
            "MOPAC",
            "NET ATOMIC CHARGES",
            "FINAL HEAT OF FORMATION",
            "CARTESIAN COORDINATES",
            "Empirical Formula"
        ],

        "orca": [
            "O   R   C   A",
            "Program Version",
            "MULLIKEN ATOMIC CHARGES",
            "LOEWDIN ATOMIC CHARGES",
            "CHELPG Charges",
            "TOTAL SCF ENERGY"
        ]
    }

    scores = {
        "xtb": 0,
        "mopac": 0,
        "orca": 0
    }

    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    for software, patterns in signatures.items():
        for pattern in patterns:
            if pattern in content:
                scores[software] += 1

    best_match = max(scores, key=scores.get)

    if scores[best_match] == 0:
        return "unknown"

    return best_match

def chrg_file_parser (_file = None, _type = None):
    """ Function doc 
    
    for index, chg in enumerate(input_files['charges']):
        system.mmState.charges[index] = chg
    
    """
    _type = _file.split('.')
    
    print('_type', _type)
    if len(_type) > 1:
        _type = _type[-1]
    else:
        _type = 'unk'
    
    
    if _type == 'MOL2' or _type == 'mol2':
        charges = extract_mol2_charges(_file)
    
    
    elif _type =='out' or _type == 'log':
        
        out_file_type = detect_qm_log_type(_file)
        
        if out_file_type == 'orca':
            charges = extract_orca_charges(_file)
            #pass
        elif out_file_type == 'xtb':
            charges = extract_xtb_charges(_file)
        
        elif out_file_type == 'mopac':
            charges = extract_mopac_charges(_file)


    elif _type == 'txt' or _type == 'unk':
        charges = extract_txt_charges (_file)
        
    elif _type == 'aux':
        charges = extract_aux_charges(_file)
    
    else:
        pass
    
    #print (charges)
    return charges







#chrg_file_parser (_file = '/home/fernando/programs/pDynamo3/scratch/bala/5-molsys_OPLS_ORCA QC MODEL_QC22_single_point.orca.log')
#chrg_file_parser (_file = '/home/fernando/programs/pDynamo3/scratch/MOPACJob.out')
#chrg_file_parser (_file = '/home/fernando/programs/pDynamo3/scratch/XTBJob.log')
#chrg_file_parser (_file = '/home/fernando/programs/pDynamo3/scratch/5-molsys_OPLS_ORCA QC MODEL_QC22_single_point.orca.log')

'''
files =[
'/home/fernando/programs/pDynamo3/scratch/4-molsys_OPLS_ORCA QC MODEL_QC22_single_point.orca.log',
'/home/fernando/programs/pDynamo3/scratch/XTBJob.log',
#'/home/fernando/programs/pDynamo3/scratch/XTBJob.log',
'/home/fernando/programs/pDynamo3/scratch/MOPACJob.out',
]

for _file in files:
    print(detect_qm_log_type(_file))
    



data = extract_orca_charges (filename = '/home/fernando/programs/pDynamo3/scratch/4-molsys_OPLS_ORCA QC MODEL_QC22_single_point.orca.log')
data = extract_orca_charges (filename = '/home/fernando/programs/pDynamo3/scratch/5-molsys_OPLS_ORCA QC MODEL_QC22_single_point.orca.log')
#data = read_XTB_log  (filename ='/home/fernando/programs/pDynamo3/scratch/XTBJob.log')
#data = extract_xtb_charges('/home/fernando/programs/pDynamo3/scratch/XTBJob.log')
#data = extract_mopac_charges('/home/fernando/programs/pDynamo3/scratch/MOPACJob.out')
data = extract_aux_charges('/home/fernando/programs/Gabedit64/maria01.aux')

print (data)
#read_MOL2 (filein = "/home/fernando/Desktop/NoName.mol2")
'''
