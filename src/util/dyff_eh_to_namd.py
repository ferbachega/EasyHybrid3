#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: DYFF / DYFF-1.0-eh system (.pkl or .yaml) -> CHARMM-format PSF + PRM +
#          PDB (+ CRD) and a NAMD configuration, with an energy check against pDynamo.
#
#  Copyright 2022-2026 Fernando Bachega
#  GNU General Public License v3 or later (see the other EasyHybrid modules).
#
#  2026-10-15. Successor of dyff_to_charmm.py for systems built with EasyHybrid's
#  DYFF (residue library + xTB-derived terms: every bond and angle has its own
#  parameters, torsions are general cosine series, polar hydrogens and waters have
#  their own LJ types). What changes with respect to dyff_to_charmm.py:
#
#    ATOM TYPES. CHARMM looks parameters up by atom-type tuple. Here every atom
#    starts from its LJ type and the classes are refined (as a Weisfeiler-Lehman
#    colouring) by the parameters of the terms each atom takes part in and the
#    classes of its partners in those terms, until stable; any tuple still mapped
#    to two parameter sets splits its atoms. Result: the smallest set of types for
#    which every bond/angle/dihedral/improper type tuple has ONE parameter set
#    (equal residues share types; a ligand gets its own).
#
#    TORSIONS, EXACT. sum_n c_n cos(n phi) -> one CHARMM line K (1 + cos(n phi -
#    delta)) per n (K = |c_n|, delta = 0 / 180 deg; NAMD takes several lines for
#    the same four types). Only a constant per dihedral is lost (reported).
#
#    LENNARD-JONES, EXACT. DYFF combines sigma AND epsilon geometrically (OPLS
#    form); NAMD does the same with `vdwGeometricSigma yes` (written in the
#    configuration), so no NBFIX table is needed. 1-4 pairs: pDynamo's own 1-4
#    LJ parameters (already scaled) in the 1-4 columns of the PRM, electrostatics
#    scaled by electrostaticScale14 (NAMD `1-4scaling`).
#
#    ANGLES AND OUT-OF-PLANE: APPROXIMATE. DYFF bends and out-of-plane terms are
#    cosine series; CHARMM/NAMD only has harmonic angles and harmonic impropers.
#    As in dyff_to_charmm.py: the harmonic with the same minimum and the same
#    curvature (d2E/dx2) -- equal near equilibrium, different for large
#    deformations; the out-of-plane angle omega is turned into the improper
#    dihedral psi through the geometric ratio d omega / d psi (mean over the
#    instances of a type tuple). The size of the approximation is reported by
#    validate().
#
#  Usage:
#      python3 dyff_eh_to_namd.py system.pkl output/prefix [--namd /path/namd3] [--no-validate]
#
import math
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict

try:
    from util.dyff_to_charmm import (kj_to_kcal, find_minimum, oop_type_curvature, geometric_domega_dpsi,
                                     charmm_psi, dyff_omega, _series_value_and_derivatives)
except ImportError:                                            # . run as a script from src/util
    from dyff_to_charmm import (kj_to_kcal, find_minimum, oop_type_curvature, geometric_domega_dpsi,
                                charmm_psi, dyff_omega, _series_value_and_derivatives)

_LJ_2_1_6 = 2.0 ** (1.0 / 6.0)
_ROUND = 8                                                     # . digits used to compare parameters


#=====================================================================================
# . Loading.
#=====================================================================================
def load_system(path):
    """ A pDynamo System from a .pkl or .yaml (EasyHybrid's src/ must be importable:
        the DYFF MM model class of EasyHybrid lives in pdynamo/dyff_edits.py). """
    import pBabel  # noqa: F401 -- import order: avoids a circular import in pMolecule
    from pBabel import ImportSystem
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        return ImportSystem(path, log=None)


def parse_atom_path(path):
    """ pDynamo atom path "SEGMENT:RESNAME.NUMBER:ATOM" (the segment may be empty)
        -> ( segment, residue name, residue number, atom name ). """
    parts = (path or "").split(":")
    atom = parts[-1] if parts and parts[-1] else "X"
    residue = parts[-2] if len(parts) >= 2 else "UNK.1"
    segment = parts[-3] if len(parts) >= 3 and parts[-3] else "SYS"
    name, _, number = residue.partition(".")
    m = re.match(r"-?\d+", number or "")
    return (segment, name or "UNK", int(m.group(0)) if m else 1, atom)


#=====================================================================================
# . Geometry helpers.
#=====================================================================================
def _angle(a, b, c):
    u = [a[k] - b[k] for k in range(3)]
    v = [c[k] - b[k] for k in range(3)]
    nu = math.sqrt(sum(x * x for x in u)); nv = math.sqrt(sum(x * x for x in v))
    cos = max(-1.0, min(1.0, sum(u[k] * v[k] for k in range(3)) / (nu * nv)))
    return math.acos(cos)


def _wrap(x):
    """ angle difference into (-pi, pi] """
    while x > math.pi: x -= 2.0 * math.pi
    while x <= -math.pi: x += 2.0 * math.pi
    return x


#=====================================================================================
# . Extraction.
#=====================================================================================
class NAMDData:
    pass


def _type_name(k):
    """ 'D' + base 36 (at most 6 characters in total) """
    digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    s = ""
    while True:
        s = digits[k % 36] + s
        k //= 36
        if k == 0: break
    return "D" + s


def _canonical(classes, palindromic=True):
    t = tuple(classes)
    return min(t, t[::-1]) if palindromic else t


def extract(system):
    """ Everything the writers need (NAMDData), from a system with a DYFF MM model. """
    from pScientific import PeriodicTable
    ms = system.mmState
    atoms = system.atoms
    n = len(atoms)
    xyz = [(system.coordinates3[i][0], system.coordinates3[i][1], system.coordinates3[i][2]) for i in range(n)]
    data = NAMDData()
    data.n = n
    data.coords = xyz
    data.residues = [parse_atom_path(atoms[i].path) for i in range(n)]
    data.elementSymbols = [PeriodicTable[atoms[i].atomicNumber].symbol for i in range(n)]
    data.masses = [PeriodicTable[atoms[i].atomicNumber].mass for i in range(n)]
    data.charges = [float(q) for q in ms.charges]
    data.dyffTypes = [ms.atomTypes[k] for k in ms.atomTypeIndices]
    data.electrostaticScale14 = float(system.mmModel.electrostaticScale14)
    # . 2026-10-15: containers by class (and label: pDynamo's OPLS has two
    #   FourierDihedralContainers, "Fourier Dihedral" and "Fourier Out-Of-Plane", and
    #   a second HarmonicBondContainer, "Urey-Bradley", the H-H bond of TIP3P water)
    supported = {"HarmonicBondContainer", "HarmonicAngleContainer", "CosineAngleContainer", "CosineDihedralContainer",
                 "FourierDihedralContainer", "CosineOutOfPlaneContainer"}
    unknown = sorted({type(c).__name__ for c in ms.mmTerms} - supported)
    if unknown:
        raise ValueError("MM terms that cannot be written in CHARMM format: {}".format(", ".join(unknown)))
    by_class = defaultdict(list)
    for c in ms.mmTerms:
        by_class[type(c).__name__].append((getattr(c, "label", "") or "", c.__getstate__()))

    # . terms: (family, atoms, parameter key); families B / A / D / I (harmonic improper) / P (periodic improper)
    terms = []
    data.constant_kcal = 0.0                                   # . pDynamo energy - NAMD energy, exact parts
    # . bonds (and Urey-Bradley H-H): E = k (r - r0)^2 in both programs
    for label, st in by_class["HarmonicBondContainer"]:
        params = [(float(eq), kj_to_kcal(float(fc))) for (eq, fc) in st["parameters"]]
        for i, j, t, active in st["terms"]:
            if active:
                b0, kb = params[t]
                terms.append(("B", (i, j), (round(b0, _ROUND), round(kb, _ROUND))))
    # . harmonic angles (OPLS): E = k (theta - theta0)^2, exact
    data.angleSeries = {}
    for label, st in by_class["HarmonicAngleContainer"]:
        params = [(math.degrees(float(t0)), kj_to_kcal(float(k))) for (t0, k) in st["parameters"]]
        for i, j, k, t, active in st["terms"]:
            if active:
                th0, kth = params[t]
                terms.append(("A", (i, j, k), (round(th0, 6), round(kth, 6))))
    # . DYFF cosine angles: curvature-matched harmonic (approximation)
    for label, st in by_class["CosineAngleContainer"]:
        harm = []
        for series in st["parameters"]:
            series = [tuple(p) for p in series]
            theta0, d2e = find_minimum(series)
            harm.append((math.degrees(theta0), kj_to_kcal(d2e / 2.0), series))
        for indices, t, active in st["terms"]:
            if active:
                th0, kth, series = harm[t]
                key = (round(th0, 6), round(kth, 6))
                terms.append(("A", tuple(indices), key))
                data.angleSeries[tuple(indices)] = series
    # . dihedrals, exact: per quadruple, lines K (1 + cos(n phi - delta)) summed by (n, delta).
    #   DYFF cosine series sum_n c_n cos(n phi) -> K = |c_n|, delta 0/180 (+ a constant);
    #   OPLS Fourier terms are already K (1 + cos(n phi - delta)).
    summed = defaultdict(lambda: defaultdict(float))
    def canonical_quad(q):
        q = tuple(q)
        return q[::-1] if q[::-1] < q else q
    for label, st in by_class["CosineDihedralContainer"]:
        for indices, t, active in st["terms"]:
            if not active: continue
            q = canonical_quad(indices)
            for c, p in st["parameters"][t]:
                c, p = kj_to_kcal(float(c)), int(p)
                if p == 0:
                    data.constant_kcal += c                       # . constant: lost in NAMD
                elif abs(c) > 1.0e-12:
                    summed[q][(p, 0.0 if c > 0.0 else 180.0)] += abs(c)
                    data.constant_kcal -= abs(c)                  # . K (1 + cos(...)) adds +K
    periodic_improper = defaultdict(lambda: defaultdict(float))
    for label, st in by_class["FourierDihedralContainer"]:
        improper = "out-of-plane" in label.lower() or "improper" in label.lower()
        for i, j, k, l, t, active in st["terms"]:
            if not active: continue
            K, period, delta = st["parameters"][t]
            key = (int(period), round(math.degrees(float(delta)) % 360.0, 6))
            if improper:
                periodic_improper[(i, j, k, l)][key] += kj_to_kcal(float(K))
            else:
                summed[canonical_quad((i, j, k, l))][key] += kj_to_kcal(float(K))
    for q, series in summed.items():
        lines = tuple((round(K, _ROUND), per, d) for (per, d), K in sorted(series.items()) if abs(K) > 1.0e-12)
        if lines:
            terms.append(("D", q, lines))
    for q, series in periodic_improper.items():
        lines = tuple((round(K, _ROUND), per, d) for (per, d), K in sorted(series.items()) if abs(K) > 1.0e-12)
        if lines:
            terms.append(("P", q, lines))
    # . out-of-plane -> harmonic impropers (approximation); psi0 0 or 180 by the geometry
    data.oop = []
    for label, st in by_class["CosineOutOfPlaneContainer"]:
        curv = [oop_type_curvature([tuple(p) for p in series]) for series in st["parameters"]]
        for indices, t, active in st["terms"]:
            if not active: continue
            i, j, k, l = indices
            psi = charmm_psi(xyz[i], xyz[j], xyz[k], xyz[l])
            ratio = geometric_domega_dpsi(xyz[i], xyz[j], xyz[k], xyz[l])
            omega0 = curv[t][0]
            if abs(omega0 - 0.5 * math.pi) < 1.0e-3 or ratio is None:
                # . planar centre: the improper minimum at psi = 0 or 180
                psi0 = 0.0 if math.cos(psi) > 0.0 else 180.0
            else:
                # . pyramidal centre (e.g. amine N:Tet, omega0 != 90 deg): the DYFF minima are
                #   omega0 and pi - omega0; the nearest one to the present geometry, mapped to
                #   psi with this instance's d omega / d psi (keeps the chirality, the
                #   harmonic improper has one well only)
                omega = dyff_omega(xyz[i], xyz[j], xyz[k], xyz[l])
                target = min((omega0, math.pi - omega0), key=lambda w: abs(w - omega))
                psi0 = math.degrees(_wrap(psi + (target - omega) / ratio))
                psi0 = float(round(psi0))
            terms.append(("I", tuple(indices), (t, psi0)))
            data.oop.append((tuple(indices), t, ratio, [tuple(p) for p in st["parameters"][t]], curv[t][1]))

    # . atom classes: LJ type + the parameters of the atom's own terms (family,
    #   position in the term, parameters) -- a LOCAL signature: equal residues get
    #   equal types. (Refining also by the partners' classes, Weisfeiler-Lehman
    #   style, would end up telling every atom of a chain apart by its distance
    #   to the chain ends.) Any type tuple still mapped to two parameter sets
    #   splits the atoms of the terms involved, refined by their partners.
    lj_index = [int(k) for k in ms.ljTypeIndices]
    by_atom = defaultdict(list)
    for m, (fam, ats, key) in enumerate(terms):
        for a in ats:
            by_atom[a].append(m)
    def local(a):
        items = []
        for m in by_atom[a]:
            fam, ats, key = terms[m]
            pos = ats.index(a)
            if fam not in ("I", "P"):
                pos = min(pos, len(ats) - 1 - pos)            # . palindromic families: direction-free
            items.append((fam, pos, key))
        return ("L", lj_index[a], tuple(sorted(items, key=repr)))
    cls = [local(a) for a in range(n)]
    split_round = 0
    while True:
        seen = {}
        conflict_atoms = set()
        for fam, ats, key in terms:
            tup = (fam, _canonical([cls[b] for b in ats], fam not in ("I", "P")))
            if tup in seen and seen[tup] != key:
                conflict_atoms.update(ats)
            seen.setdefault(tup, key)
        if not conflict_atoms:
            break
        split_round += 1
        if split_round > 20:
            cls = [("U", a) if a in conflict_atoms else cls[a] for a in range(n)]   # . last resort: unique
            continue
        new = list(cls)
        for a in conflict_atoms:
            partners = tuple(sorted((terms[m][0], tuple(cls[b] for b in terms[m][1])) for m in by_atom[a]))
            new[a] = ("S", cls[a], partners)
        cls = new
    order = {}
    for c in cls:
        order.setdefault(c, len(order))
    data.typeOfAtom = [_type_name(order[c]) for c in cls]
    data.numberOfTypes = len(order)
    data.split_rounds = split_round

    # . nonbonded parameters per type (from the atom's LJ type)
    lj = ms.ljParameters.__getstate__()
    lj14 = ms.ljParameters14.__getstate__()
    data.nonbonded = {}
    data.typeInfo = {}
    for a in range(n):
        name = data.typeOfAtom[a]
        k = lj_index[a]
        data.nonbonded.setdefault(name, (kj_to_kcal(lj["epsilons"][k]), _LJ_2_1_6 * lj["sigmas"][k] / 2.0,
                                         kj_to_kcal(lj14["epsilons"][k]), _LJ_2_1_6 * lj14["sigmas"][k] / 2.0))
        seg, res, num, atom = data.residues[a]
        data.typeInfo.setdefault(name, "{} {} DYFF {} LJ {}".format(res, atom, data.dyffTypes[a], lj["parameterKeys"][k]))

    # . topology lists and parameter tables
    data.bonds, data.angles, data.dihedrals, data.impropers = [], [], [], []
    data.bondParams, data.angleParams, data.dihedralParams, data.improperParams = {}, {}, {}, {}
    data.improperPeriodic = {}                                 # . OPLS Fourier impropers, exact
    oop_ratios = defaultdict(list)
    oop_info = {q: (t, ratio, curvature) for (q, t, ratio, series, curvature) in data.oop}
    for fam, ats, key in terms:
        names = tuple(data.typeOfAtom[b] for b in ats)
        if fam == "B":
            data.bonds.append(ats); data.bondParams.setdefault(_canonical(names), key)
        elif fam == "A":
            data.angles.append(ats); data.angleParams.setdefault(_canonical(names), key)
        elif fam == "D":
            data.dihedrals.append(ats)
            tup = _canonical(names)
            if tup != names:                                   # . lines are written for the canonical order
                pass
            data.dihedralParams.setdefault(tup, key)
        elif fam == "P":
            data.impropers.append(ats); data.improperPeriodic.setdefault(names, key)
        else:
            data.impropers.append(ats)
            t, ratio, curvature = oop_info[ats]
            if ratio is not None:
                oop_ratios[names].append(ratio)
            data.improperParams.setdefault(names, (t, key[1], curvature))
    data.improperK = {}
    for names, (t, psi0, curvature) in data.improperParams.items():
        r = oop_ratios.get(names)
        mean = sum(r) / len(r) if r else 0.0
        data.improperK[names] = (kj_to_kcal(curvature * mean * mean / 2.0), psi0)

    data.box = None
    if system.symmetryParameters is not None:
        sp = system.symmetryParameters
        data.box = (sp.a, sp.b, sp.c, sp.alpha, sp.beta, sp.gamma)
    return data


#=====================================================================================
# . Writers.
#=====================================================================================
def _int_lines(values, per_line):
    line = []
    for v in values:
        line.append("{:10d}".format(v))
        if len(line) == per_line:
            yield "".join(line); line = []
    if line:
        yield "".join(line)


def write_psf(path, data, title="DYFF system exported by EasyHybrid (dyff_eh_to_namd.py)"):
    """ XPLOR-format PSF (atom types as names), EXT field widths. """
    L = ["PSF EXT", "", "{:10d} !NTITLE".format(1), "* " + title, ""]
    L.append("{:10d} !NATOM".format(data.n))
    for i in range(data.n):
        seg, res, num, atom = data.residues[i]
        L.append("{:10d} {:<8s} {:<8d} {:<8s} {:<8s} {:<6s} {:14.6f} {:13.4f} {:11d}".format(
            i + 1, seg, num, res[:8], atom[:8], data.typeOfAtom[i], data.charges[i], data.masses[i], 0))
    L.append("")
    for label, items, per in (("!NBOND: bonds", data.bonds, 4), ("!NTHETA: angles", data.angles, 3),
                              ("!NPHI: dihedrals", data.dihedrals, 2), ("!NIMPHI: impropers", data.impropers, 2)):
        L.append("{:10d} {}".format(len(items), label))
        L.extend(_int_lines([a + 1 for t in items for a in t], per * len(items[0]) if items else 8))
        L.append("")
    L += ["{:10d} !NDON: donors".format(0), "", "{:10d} !NACC: acceptors".format(0), "",
          "{:10d} !NNB".format(0), ""]
    L.extend(_int_lines([0] * data.n, 8))
    L += ["", "{:10d}{:10d} !NGRP NST2".format(1, 0), "{:10d}{:10d}{:10d}".format(0, 0, 0), "",
          "{:10d}{:10d} !NUMLP NUMLPH".format(0, 0), ""]
    open(path, "w").write("\n".join(L) + "\n")


def write_prm(path, data):
    L = ["* DYFF parameters exported by EasyHybrid (dyff_eh_to_namd.py).",
         "* NAMD: vdwGeometricSigma yes; 1-4scaling {:.4f}".format(data.electrostaticScale14),
         "* Bonds, dihedrals and Lennard-Jones: exact. Angles and impropers: harmonic with the",
         "* minimum and curvature of the DYFF cosine terms (approximation).",
         "* LJ combining: geometric sigma and epsilon -> NAMD needs 'vdwGeometricSigma yes'.",
         "*", "",
         "! atom types: name  residue atom  DYFF type  LJ type"]
    for name in sorted(data.typeInfo, key=lambda s: (len(s), s)):
        L.append("! {:<7s} {}".format(name, data.typeInfo[name]))
    L += ["", "BONDS"]
    for (a, b), (b0, kb) in data.bondParams.items():
        L.append("{:<6s} {:<6s} {:14.6f} {:12.6f}".format(a, b, kb, b0))
    L += ["", "ANGLES"]
    for (a, b, c), (th0, kth) in data.angleParams.items():
        L.append("{:<6s} {:<6s} {:<6s} {:14.6f} {:12.6f}".format(a, b, c, kth, th0))
    L += ["", "DIHEDRALS"]
    for (a, b, c, d), lines in data.dihedralParams.items():
        for K, nper, delta in lines:                           # . consecutive lines = multiple terms
            L.append("{:<6s} {:<6s} {:<6s} {:<6s} {:14.8f} {:3d} {:8.2f}".format(a, b, c, d, K, nper, delta))
    L += ["", "IMPROPERS"]
    for (a, b, c, d), (K, psi0) in data.improperK.items():
        # . NAMD measures the improper dihedral with the opposite sign of charmm_psi()
        #   (checked by a single point: only matters for psi0 other than 0 / 180)
        namd_psi0 = -psi0 if abs(psi0) not in (0.0, 180.0) else psi0
        L.append("{:<6s} {:<6s} {:<6s} {:<6s} {:14.6f} {:3d} {:8.2f}".format(a, b, c, d, K, 0, namd_psi0))
    for (a, b, c, d), lines in data.improperPeriodic.items():
        # . periodic impropers K (1 + cos(n psi - delta)); same sign rule as above
        for K, nper, delta in lines:
            namd_delta = (-delta) % 360.0 if delta % 180.0 else delta
            L.append("{:<6s} {:<6s} {:<6s} {:<6s} {:14.8f} {:3d} {:8.2f}".format(a, b, c, d, K, nper, namd_delta))
    L += ["", "NONBONDED nbxmod 5 atom cdiel shift vatom vdistance vswitch -",
          "cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5",
          "!type    ignored   epsilon      Rmin/2    ignored  eps,1-4   Rmin/2,1-4"]
    for name in sorted(data.nonbonded, key=lambda s: (len(s), s)):
        eps, rh, eps14, rh14 = data.nonbonded[name]
        L.append("{:<6s} {:8.4f} {:12.8f} {:12.8f} {:8.4f} {:12.8f} {:12.8f}".format(name, 0.0, -eps, rh, 0.0, -eps14, rh14))
    L += ["", "END", ""]
    open(path, "w").write("\n".join(L))


_NAMD_LINE = re.compile(r"^\*\s*NAMD:\s*vdwGeometricSigma yes;\s*1-4scaling\s+([0-9.]+)")


def namd_settings_from_prm(paths):
    """ {"vdwGeometricSigma": True, "1-4scaling": x} when one of the CHARMM parameter
        files was written by this module (its title declares them), else {} -- used
        by namd_runner.build_namd_config so the "Run NAMD" window gets them right. """
    if not paths:
        return {}
    if isinstance(paths, str):
        paths = [paths]
    for path in paths:
        try:
            with open(path) as handle:
                for _, line in zip(range(5), handle):
                    m = _NAMD_LINE.match(line)
                    if m:
                        return {"vdwGeometricSigma": True, "1-4scaling": float(m.group(1))}
        except OSError:
            continue
    return {}


def write_pdb(path, data):
    L = []
    for i in range(data.n):
        seg, res, num, atom = data.residues[i]
        x, y, z = data.coords[i]
        name = atom[:4] if len(atom) >= 4 else " " + atom          # . PDB convention for 1-3 letter names
        L.append("ATOM  {:>5d} {:<4s} {:>3s} {:1s}{:>4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}{:>2s}".format(
            (i + 1) % 100000, name, res[:3], seg[:1], num % 10000, x, y, z, 1.0, 0.0, seg[:4], data.elementSymbols[i]))
    L.append("END")
    open(path, "w").write("\n".join(L) + "\n")


def write_crd(path, data, title="DYFF system exported by EasyHybrid (dyff_eh_to_namd.py)"):
    """ CHARMM extended coordinate file (for CHARMM itself; NAMD reads the PDB). """
    L = ["* " + title, "*", "{:10d}  EXT".format(data.n)]
    for i in range(data.n):
        seg, res, num, atom = data.residues[i]
        x, y, z = data.coords[i]
        L.append("{:10d}{:10d}  {:<8s}  {:<8s}{:20.10f}{:20.10f}{:20.10f}  {:<8s}  {:<8s}{:20.10f}".format(
            i + 1, num, res[:8], atom[:8], x, y, z, seg[:8], str(num), 0.0))
    open(path, "w").write("\n".join(L) + "\n")


def write_namd_config(path, prefix, data, run="singlepoint", cutoff=None):
    """ NAMD input. run = "singlepoint" (run 0, no cutoff: energy check) or
        "minimize" (1000 steps). The DYFF-specific lines (vdwGeometricSigma,
        1-4scaling) must be kept in any other NAMD input made for these files. """
    base = os.path.basename(prefix)
    if cutoff is None:
        cutoff = 12.0 if data.box else 999.0
    L = ["# NAMD input for {} (EasyHybrid DYFF -> CHARMM format, dyff_eh_to_namd.py)".format(base),
         "structure          {}.psf".format(base),
         "coordinates        {}.pdb".format(base),
         "paraTypeCharmm     on",
         "parameters         {}.prm".format(base),
         "",
         "# . DYFF specific: geometric combining of sigma (and epsilon), 1-4 electrostatics",
         "vdwGeometricSigma  yes",
         "exclude            scaled1-4",
         "1-4scaling         {:.4f}   ;# NAMD 3 also accepts the newer name oneFourScaling".format(data.electrostaticScale14),
         "",
         "temperature        0" if run == "singlepoint" else "temperature        300",
         "cutoff             {:.1f}".format(cutoff),
         "pairlistdist       {:.1f}".format(cutoff + 2.0),
         "switching          off" if not data.box else "switching          on\nswitchdist         10.0"]
    if data.box:
        a, b, c = data.box[:3]
        L += ["cellBasisVector1   {:.4f} 0 0".format(a), "cellBasisVector2   0 {:.4f} 0".format(b),
              "cellBasisVector3   0 0 {:.4f}".format(c), "wrapAll            on", "PME                yes",
              "PMEGridSpacing     1.0"]
    L += ["", "outputName         {}_{}".format(base, run), "binaryOutput       no",
          "outputEnergies     {}".format(1 if run == "singlepoint" else 100), ""]
    L.append("run 0" if run == "singlepoint" else "minimize 1000")
    open(path, "w").write("\n".join(L) + "\n")


#=====================================================================================
# . Energy check.
#=====================================================================================
def pdynamo_energies(system):
    """ kcal/mol per CHARMM component. Non periodic: no cutoff (NBModelFull), every
        component. Periodic: the system's own NB model, and only the bonded
        components are meaningful for a comparison (cutoff / PME differ). """
    import contextlib
    import io
    from pMolecule.NBModel import NBModelFull
    s = system
    old = s.nbModel
    periodic = s.symmetryParameters is not None
    if not periodic:
        s.DefineNBModel(NBModelFull.WithDefaults())
    with contextlib.redirect_stdout(io.StringIO()):
        s.Energy(log=None)
    t = {k: kj_to_kcal(v) for k, v in s.scratch.energyTerms.items()}
    if old is not None and not periodic:
        s.DefineNBModel(old)
    g = lambda *keys: sum(t.get(k, 0.0) for k in keys)
    return {"BOND": g("Harmonic Bond", "Urey-Bradley"), "ANGLE": g("Cosine Angle", "Harmonic Angle"),
            "DIHED": g("Cosine Dihedral", "Fourier Dihedral"), "IMPRP": g("Cosine Out-Of-Plane", "Fourier Out-Of-Plane"),
            "ELECT": g("MM/MM Electrostatic", "MM/MM 1-4 Electrostatic"),
            "VDW": g("MM/MM Lennard-Jones", "MM/MM 1-4 Lennard-Jones")}


def harmonic_estimates(data):
    """ What the approximate terms should give in NAMD (kcal/mol), evaluated here
        from the written parameters -- separates "file written correctly" (NAMD ==
        this) from "approximation" (this vs pDynamo, after removing constants). """
    x = data.coords
    ang = 0.0
    ang_dyff_shift = 0.0
    for ats in data.angles:
        th0, kth = data.angleParams[_canonical(tuple(data.typeOfAtom[b] for b in ats))]
        th = _angle(x[ats[0]], x[ats[1]], x[ats[2]])
        ang += kth * (th - math.radians(th0)) ** 2
        series = data.angleSeries.get(ats)
        if series is not None:                                  # . DYFF cosine angle (harmonic angles: exact)
            e0, _, _ = _series_value_and_derivatives(series, math.radians(th0))
            ang_dyff_shift += kj_to_kcal(e0)                    # . DYFF energy at its own minimum
    imp = 0.0
    for ats in data.impropers:
        names = tuple(data.typeOfAtom[b] for b in ats)
        psi = charmm_psi(*[x[b] for b in ats])
        if names in data.improperK:
            K, psi0 = data.improperK[names]
            imp += K * _wrap(psi - math.radians(psi0)) ** 2
        else:
            imp += sum(K * (1.0 + math.cos(n * psi - math.radians(d))) for K, n, d in data.improperPeriodic[names])
    return {"ANGLE": ang, "IMPRP": imp, "angle_minimum_energy": ang_dyff_shift}


def find_namd():
    for c in (os.getenv("NAMD"), shutil.which("namd3"), shutil.which("namd2")):
        if c and os.access(c, os.X_OK):
            return c
    import glob
    for c in sorted(glob.glob(os.path.expanduser("~/programs/NAMD*/namd3")) +
                    glob.glob(os.path.expanduser("~/programs/NAMD*/namd2")), reverse=True):
        if os.access(c, os.X_OK):
            return c
    return None


def run_namd(config, namd=None, cpus=1):
    """ ENERGY components of the first ENERGY: line (kcal/mol). """
    namd = namd or find_namd()
    if namd is None:
        raise RuntimeError("NAMD not found (set $NAMD or put namd3/namd2 in the PATH)")
    folder, name = os.path.split(os.path.abspath(config))
    result = subprocess.run([namd, "+p{}".format(cpus), name], cwd=folder, capture_output=True, text=True, timeout=3600)
    log = result.stdout + result.stderr
    open(os.path.splitext(config)[0] + ".log", "w").write(log)
    titles = None
    for line in log.splitlines():
        if line.startswith("ETITLE:"):
            titles = line.split()[1:]
        elif line.startswith("ENERGY:") and titles:
            values = line.split()[1:]
            return {k: float(v) for k, v in zip(titles, values)}
    err = [l for l in log.splitlines() if "ERROR" in l]
    raise RuntimeError("NAMD failed: " + ("; ".join(err[:3]) if err else log[-800:]))


def validate(system, data, prefix, namd=None):
    """ Single point in NAMD (no cutoff) vs pDynamo (NBModelFull). Returns rows
        (component, pDynamo, NAMD, expected NAMD, difference, comment). """
    config = prefix + "_singlepoint.namd"
    write_namd_config(config, prefix, data, run="singlepoint")
    namd_e = run_namd(config, namd)
    # . both programs on the same coordinates: the PDB has 3 decimals
    import copy
    s = copy.deepcopy(system)
    for i in range(data.n):
        for k in range(3):
            s.coordinates3[i, k] = round(float(s.coordinates3[i, k]), 3)
    pd = pdynamo_energies(s)
    rounded = copy.copy(data)
    rounded.coords = [tuple(round(v, 3) for v in c) for c in data.coords]
    est = harmonic_estimates(rounded)
    rows = []
    rows.append(("BOND", pd["BOND"], namd_e["BOND"], pd["BOND"], "exact"))
    rows.append(("DIHED", pd["DIHED"], namd_e["DIHED"], pd["DIHED"] - data.constant_kcal,
                 "exact (pDynamo minus the constants of the cosine series, {:+.4f})".format(data.constant_kcal)))
    if not data.box:                                        # . periodic: PME / cutoff differ, bonded terms only
        rows.append(("ELECT", pd["ELECT"], namd_e["ELECT"], pd["ELECT"], "exact (incl. 1-4)"))
        rows.append(("VDW", pd["VDW"], namd_e["VDW"], pd["VDW"], "exact (incl. 1-4)"))
    rows.append(("ANGLE", pd["ANGLE"], namd_e["ANGLE"], est["ANGLE"],
                 "exact (harmonic)" if not data.angleSeries else
                 "harmonic approximation; pDynamo relative to the minima: {:.4f}".format(pd["ANGLE"] - est["angle_minimum_energy"])))
    rows.append(("IMPRP", pd["IMPRP"], namd_e["IMPRP"], est["IMPRP"],
                 "exact (Fourier)" if not data.improperK else "harmonic approximation"))
    return rows


#=====================================================================================
# . Orchestration.
#=====================================================================================
def convert(path_or_system, prefix, crd=True):
    system = load_system(path_or_system) if isinstance(path_or_system, str) else path_or_system
    data = extract(system)
    folder = os.path.dirname(os.path.abspath(prefix))
    os.makedirs(folder, exist_ok=True)
    write_psf(prefix + ".psf", data)
    write_prm(prefix + ".prm", data)
    write_pdb(prefix + ".pdb", data)
    if crd:
        write_crd(prefix + ".crd", data)
    write_namd_config(prefix + "_minimize.namd", prefix, data, run="minimize")
    return system, data


def summary(data):
    return {"atoms": data.n, "atom_types": data.numberOfTypes, "bonds": len(data.bonds), "angles": len(data.angles),
            "dihedrals": len(data.dihedrals), "impropers": len(data.impropers),
            "bond_types": len(data.bondParams), "angle_types": len(data.angleParams),
            "dihedral_types": len(data.dihedralParams), "improper_types": len(data.improperK) + len(data.improperPeriodic),
            "electrostaticScale14": data.electrostaticScale14, "box": data.box}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="DYFF (EasyHybrid) system -> NAMD (CHARMM-format PSF/PRM/PDB/CRD)")
    parser.add_argument("system", help=".pkl or .yaml")
    parser.add_argument("prefix", help="output prefix (folder/name)")
    parser.add_argument("--namd", default=None)
    parser.add_argument("--no-validate", action="store_true")
    args = parser.parse_args()
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if here not in sys.path: sys.path.insert(0, here)
    system, data = convert(args.system, args.prefix)
    for k, v in summary(data).items():
        print("{:>16s}: {}".format(k, v))
    if not args.no_validate:
        for comp, pd, nm, exp, note in validate(system, data, args.prefix, args.namd):
            print("{:>6s}  pDynamo {:14.4f}  NAMD {:14.4f}  expected {:14.4f}  diff {:10.4f}  {}".format(comp, pd, nm, exp, nm - exp, note))
