#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: DYFF (pDynamo's UFF implementation) -> CHARMM PSF/PRM/PDB export
#
#  Copyright 2022-2026 Fernando Bachega
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
#      Converts a pDynamo3 System pickled with the DYFF MM model (pDynamo's
#      own implementation of Rappe et al.'s Universal Force Field, UFF --
#      see ~/programs/pDynamo3/parameters/forceFields/dyff/rawData/00ReadMe.txt)
#      into a CHARMM-format PSF + PRM + PDB triple that NAMD can run
#      DIRECTLY with the real UFF energetics -- not a re-parameterisation
#      into some other force field. Once written, the 3 files plug straight
#      into EasyHybrid's existing "Run NAMD" window in its CHARMM mode
#      (structure=<prefix>.psf, parameters=[<prefix>.prm],
#      coordinates=<prefix>.pdb, one_four_scaling=<see printed report>).
#
#      THE CENTRAL PROBLEM this module solves: DYFF's bond term is a plain
#      harmonic (exactly CHARMM's own bond functional form -- an EXACT,
#      lossless conversion), but its angle and out-of-plane (improper)
#      terms are general N-term Fourier-in-cosine expansions
#      (E = sum_p c_p * cos(p*angle), confirmed against the real evaluator,
#      pMolecule/MMModel/__extensions__/csource/CosineTermEnergies.c),
#      while CHARMM/NAMD's angle and improper terms are STRICTLY harmonic
#      (K*(x-x0)**2) -- there is no NAMD keyword or PRM section that can
#      represent a general Fourier bend. This module APPROXIMATES those
#      two term types with a curvature-matched harmonic (same equilibrium
#      value, same local stiffness -- i.e. same vibrational frequency
#      around the minimum -- but a different shape away from it):
#      K_harmonic = ( d2E/dx2 at the true minimum ) / 2, which is the
#      unique harmonic whose second derivative matches DYFF's at that
#      point (a Taylor-expansion identity, not a fit).
#
#      DYFF's dihedral term, by contrast, turns out to need NO
#      approximation at all: every dihedral in this system's real data is
#      exactly a 2-term series c0 + c1*cos(n*phi) with |c0| == |c1| (the
#      signature of UFF's own E = (Vphi/2)*(1 - cos(n*phi0)*cos(n*phi))
#      form whenever cos(n*phi0) = +-1, which the UFF parameterisation
#      always arranges) -- an EXACT algebraic match to CHARMM's own
#      Kchi*(1 + cos(n*phi - delta)) form (delta = 180 deg when c1 == -c0,
#      0 deg when c1 == +c0, Kchi = |c0|). Verified against this system's
#      actual persisted parameters before relying on it, not assumed from
#      the literature.
#
#      Nonbonded: DYFF uses LJForm.OPLS, i.e. GEOMETRIC combining for
#      BOTH epsilon and sigma (matches UFF's own paper exactly). NAMD's
#      built-in CHARMM combining rule is geometric epsilon but ARITHMETIC
#      Rmin -- so every cross atom-type pair is written as an explicit
#      NBFIX override (exact geometric value), rather than relying on
#      NAMD's own default combining, which would silently be wrong for
#      any hetero-atom-type pair.
#
#      The one place true per-INSTANCE geometry is needed (not just
#      per-type constants) is converting DYFF's out-of-plane angle omega
#      (angle between the i-j bond and the normal of the j,k,l plane) into
#      CHARMM's improper dihedral psi (a totally different geometric
#      definition -- the standard 4-body torsion angle of the same 4
#      atoms). See geometric_domega_dpsi() -- a numeric, per-instance
#      finite-difference ratio, averaged over every real instance of a
#      given parameter type to get one representative Kpsi per type.
#
import math
import os
import re
import sys


#=====================================================================================
# . Unit conversion.
#=====================================================================================
_KJ_PER_KCAL = 4.184
def kj_to_kcal(value):
    return value / _KJ_PER_KCAL


#=====================================================================================
# . Loading.
#=====================================================================================
def load_system(yaml_path, pdynamo3_home=None):
    """ Loads a pDynamo3 System pickled as YAML. Importing pBabel BEFORE
        pCore.YAMLUnpickle sidesteps a real circular-import bug in this
        pDynamo3 install (pMolecule.MMModel.MMModel imports pBabel, which
        imports pMolecule back before pMolecule has finished initialising,
        if pBabel is not already fully imported first -- confirmed by
        hand, not a hypothetical).
    """
    if pdynamo3_home and pdynamo3_home not in sys.path:
        sys.path.append(pdynamo3_home)
    import pBabel  # noqa: F401 -- import order matters, see docstring.
    from pCore import YAMLUnpickle
    return YAMLUnpickle(yaml_path)


#=====================================================================================
# . Cosine-series helpers (angles / out-of-plane terms).
#=====================================================================================
def _series_value_and_derivatives(coefficientPeriodPairs, x):
    """ E, dE/dx, d2E/dx2 for E(x) = sum_p c_p * cos(p*x), x in radians. """
    e = de = d2e = 0.0
    for c, p in coefficientPeriodPairs:
        px = p * x
        e   += c * math.cos(px)
        de  += -c * p * math.sin(px)
        d2e += -c * p * p * math.cos(px)
    return e, de, d2e


def find_minimum(coefficientPeriodPairs, xLower=1.0e-4, xUpper=None, nGrid=3600):
    """ Global minimum of E(x) = sum_p c_p*cos(p*x) over (xLower, xUpper)
        (radians) by a fine grid scan (robust against multiple local
        minima/maxima, which a general N-term series can have) followed
        by a few Newton steps on dE/dx = 0 to refine it. Returns
        (xMinimum, d2E/dx2 at xMinimum).
    """
    if xUpper is None:
        xUpper = math.pi - 1.0e-4
    best_x, best_e = None, None
    for i in range(nGrid + 1):
        x = xLower + (xUpper - xLower) * i / nGrid
        e, _, _ = _series_value_and_derivatives(coefficientPeriodPairs, x)
        if best_e is None or e < best_e:
            best_e, best_x = e, x
    x = best_x
    for _ in range(50):
        _, de, d2e = _series_value_and_derivatives(coefficientPeriodPairs, x)
        if abs(d2e) < 1.0e-10:
            break
        step = de / d2e
        x -= step
        x = min(max(x, xLower), xUpper)
        if abs(step) < 1.0e-12:
            break
    _, _, d2e = _series_value_and_derivatives(coefficientPeriodPairs, x)
    return x, d2e


def angle_type_to_harmonic(coefficientPeriodPairs):
    """ (theta0Degrees, KthetaKcalPerMolPerRad2) -- see module docstring
        for why this is a curvature-matched harmonic, not an exact
        conversion.
    """
    theta0, d2e = find_minimum(coefficientPeriodPairs)
    kTheta = kj_to_kcal(d2e / 2.0)
    return math.degrees(theta0), kTheta


def dihedral_type_to_charmm(coefficientPeriodPairs):
    """ (KchiKcalPerMol, n, deltaDegrees) -- an EXACT conversion for a
        2-term (c0, c1) series with |c0| == |c1| (see module docstring);
        raises if a real dihedral parameter type doesn't fit that shape,
        rather than silently mis-converting it.
    """
    if len(coefficientPeriodPairs) != 2:
        raise ValueError("Expected a 2-term dihedral series, got {}.".format(coefficientPeriodPairs))
    (c0, p0), (c1, p1) = coefficientPeriodPairs
    if p0 != 0:
        (c0, p0), (c1, p1) = (c1, p1), (c0, p0)
    if p0 != 0 or p1 <= 0:
        raise ValueError("Unexpected dihedral series periods: {}.".format(coefficientPeriodPairs))
    if abs(abs(c0) - abs(c1)) > 1.0e-6 * max(abs(c0), 1.0):
        raise ValueError("Dihedral series is not |c0|==|c1| -- not exactly representable: {}."
                          .format(coefficientPeriodPairs))
    delta = 180.0 if (c1 < 0) == (c0 > 0) else 0.0
    # . c0 and c1 have the same magnitude; c0's sign matches the physical
    #   V/2 prefactor (always taken positive by convention -- delta
    #   already encodes the sign relationship).
    kChi = kj_to_kcal(abs(c0))
    return kChi, p1, delta


def oop_type_curvature(coefficientPeriodPairs):
    """ (omega0Radians, d2E/domega2 at omega0) for a DYFF out-of-plane
        term. Finds the true minimum generically (does not assume the
        common c0==c1, omega0=90deg shape), so it also covers the
        3-term cases DYFF uses for some heavier elements.
    """
    return find_minimum(coefficientPeriodPairs, xLower=1.0e-4, xUpper=math.pi / 2.0 - 1.0e-4)


#=====================================================================================
# . Geometry (for the omega -> psi improper conversion).
#=====================================================================================
def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

def _add(a, b, s=1.0):
    return (a[0] + b[0] * s, a[1] + b[1] * s, a[2] + b[2] * s)

def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def _norm(a):
    return math.sqrt(_dot(a, a))

def _unit(a):
    n = _norm(a)
    return (a[0] / n, a[1] / n, a[2] / n)


def dyff_omega(ri, rj, rk, rl):
    """ DYFF's own out-of-plane angle: the angle between the i-j bond and
        the normal of the (j,k,l) plane. Reproduces
        CosineTermEnergy_OutOfPlane's geometry exactly (see
        pMolecule/MMModel/__extensions__/csource/CosineTermEnergies.c).
    """
    xkj = _unit(_sub(rk, rj))
    xlj = _unit(_sub(rl, rj))
    n = _unit(_cross(xkj, xlj))
    xij = _unit(_sub(ri, rj))
    cosPhi = max(-1.0, min(1.0, _dot(n, xij)))
    return math.acos(cosPhi)


def charmm_psi(ri, rj, rk, rl):
    """ The standard 4-body torsion angle of (i,j,k,l) taken as a chain
        i-j-k-l -- exactly what CHARMM/NAMD compute for an "improper"
        (same formula as a proper dihedral; only the physical
        interpretation differs). Returned in (-pi, pi].
    """
    b1 = _sub(rj, ri)
    b2 = _sub(rk, rj)
    b3 = _sub(rl, rk)
    n1 = _cross(b1, b2)
    n2 = _cross(b2, b3)
    m1 = _cross(n1, _unit(b2))
    x = _dot(n1, n2)
    y = _dot(m1, n2)
    return math.atan2(y, x)


def geometric_domega_dpsi(ri, rj, rk, rl, step=1.0e-4):
    """ Numeric d(omega)/d(psi) at the given (real) geometry: perturbs
        atom i by +-step along the current DYFF-normal direction and
        takes the central-difference ratio of the resulting changes in
        omega and psi. This is the purely geometric conversion factor
        needed because DYFF's omega and CHARMM's psi are different
        angle definitions for the same 4 atoms (see module docstring).
    """
    xkj = _unit(_sub(rk, rj))
    xlj = _unit(_sub(rl, rj))
    n = _unit(_cross(xkj, xlj))
    ri_plus  = _add(ri, n, step)
    ri_minus = _add(ri, n, -step)
    omega_plus  = dyff_omega(ri_plus,  rj, rk, rl)
    omega_minus = dyff_omega(ri_minus, rj, rk, rl)
    psi_plus  = charmm_psi(ri_plus,  rj, rk, rl)
    psi_minus = charmm_psi(ri_minus, rj, rk, rl)
    # . Guard against the +-pi branch cut (only relevant if the real
    #   geometry sits almost exactly at psi = +-180 deg, unlikely for a
    #   near-planar sp2 improper centre but cheap to guard anyway).
    dPsi = psi_plus - psi_minus
    if dPsi > math.pi:
        dPsi -= 2.0 * math.pi
    elif dPsi < -math.pi:
        dPsi += 2.0 * math.pi
    dOmega = omega_plus - omega_minus
    if abs(dPsi) < 1.0e-10:
        return None
    return dOmega / dPsi


#=====================================================================================
# . Atom-type / label sanitising.
#=====================================================================================
def sanitize_type_label(label):
    """ DYFF atom-type labels use ':' (e.g. "C:Tet") -- replaced with '_'
        for a conservative, universally-safe CHARMM atom-type name.
    """
    return label.replace(":", "_")


_PATH_RE = re.compile(r":([^.]+)\.(\d+):(.+)")

def parse_atom_path(path):
    """ pMolecule.Atom.path is ":RESNAME.RESNUMBER:ATOMNAME" -- returns
        (resName, resNumber, atomName), or a generic fallback if a path
        is ever missing/malformed. Returned as-is (NOT truncated) --
        write_psf() separates every field with an explicit space so an
        untruncated name never risks merging with its neighbour; PDB's
        strict 4-character atom-name column truncates separately, only
        at the point write_pdb() actually needs it.
    """
    m = _PATH_RE.match(path or "")
    if m is None:
        return ("UNK", 1, "X")
    return (m.group(1), int(m.group(2)), m.group(3))


#=====================================================================================
# . Extraction.
#=====================================================================================
_LJ_2_1_6 = 2.0 ** (1.0 / 6.0)

class ChargmmData:
    """ Plain container for everything write_psf/write_prm/write_pdb need
        -- filled in by Extract().
    """
    pass


def Extract(system):
    from pScientific import PeriodicTable

    data = ChargmmData()
    atoms = system.atoms.items
    n = len(atoms)
    coords = system.coordinates3
    mmState = system.mmState

    data.n = n
    data.coords = [(coords[i][0], coords[i][1], coords[i][2]) for i in range(n)]
    data.residues = [parse_atom_path(atoms[i].path) for i in range(n)]  # (resName, resNumber, atomName)
    data.elementSymbols = [PeriodicTable[atoms[i].atomicNumber].symbol for i in range(n)]
    data.masses = [PeriodicTable[atoms[i].atomicNumber].mass for i in range(n)]
    data.charges = list(mmState.charges)

    typeTable = [sanitize_type_label(label) for label in mmState.atomTypes]
    baseTypeLabels = [typeTable[idx] for idx in mmState.atomTypeIndices]

    # . DYFF's own bond/angle parameters are BOND-ORDER dependent (e.g. a
    #   ring C:Tri-C:Tri bond and an exocyclic C:Tri-C:Tri bond to a
    #   carbonyl carbon get genuinely different lengths/force constants
    #   -- confirmed by hand against this system's real parameter table:
    #   several coarse (type,type[,type]) triples collapse two or more
    #   DIFFERENT numeric bond/angle parameter sets onto the same CHARMM
    #   type key, which is a real, silent fidelity loss since CHARMM
    #   parameter lookup only ever sees the static per-ATOM type, not the
    #   bond order). Fix: refine every atom's type with a "_D"/"_S" tag
    #   for whether ANY of its real bonds is a genuine double bond
    #   (order >= 1.9, via DYFF's own DYFFBondOrder()) -- mechanically
    #   derived from the same bond-order data DYFF itself used to pick
    #   parameters, not guessed from chemistry. Verified (see this
    #   function's own call site / development notes) to eliminate every
    #   bond and angle value collision in this system; harmless
    #   (produces an unused "_D" variant) for types that never carry a
    #   double bond, e.g. C:Tet/H/O:Tet.
    from pMolecule.MMModel.DYFFUtilities import DYFFBondOrder
    maxBondOrder = [0.0] * n
    for bond in system.connectivity.bonds:
        i, j = bond.node1.index, bond.node2.index
        bondOrder, _tag = DYFFBondOrder(i, j, system.connectivity)
        if bondOrder > maxBondOrder[i]:
            maxBondOrder[i] = bondOrder
        if bondOrder > maxBondOrder[j]:
            maxBondOrder[j] = bondOrder

    def refine(i):
        return baseTypeLabels[i] + ("_D" if maxBondOrder[i] >= 1.9 else "_S")

    data.atomTypeLabels = [refine(i) for i in range(n)]
    baseOfRefined = {}
    for i in range(n):
        baseOfRefined[data.atomTypeLabels[i]] = baseTypeLabels[i]
    refinedTypesUsed = sorted(baseOfRefined.keys())

    # . Nonbonded (LJForm.OPLS -- geometric combining for BOTH epsilon and
    #   sigma; see module docstring). eps/sigma only ever depend on the
    #   COARSE (pre-refinement) type -- both "_D"/"_S" variants of a
    #   split type share the identical physical vdW parameters, so this
    #   just looks each one up via its coarse parent.
    lj = mmState.ljParameters.__getstate__()
    lj14 = mmState.ljParameters14.__getstate__()
    ljTypeLabels = [sanitize_type_label(k) for k in lj['parameterKeys']]
    assert ljTypeLabels == typeTable, "LJ and atom-type tables are not in the same order -- unexpected."
    ljIndexOf = {label: idx for idx, label in enumerate(typeTable)}

    data.nonbondedTypes = []  # [ (label, epsKcal, rMinHalf, eps14Kcal, rMin14Half) ]
    for label in refinedTypesUsed:
        i = ljIndexOf[baseOfRefined[label]]
        epsKcal = kj_to_kcal(lj['epsilons'][i])
        rMinHalf = _LJ_2_1_6 * lj['sigmas'][i] / 2.0
        eps14Kcal = kj_to_kcal(lj14['epsilons'][i])
        rMin14Half = _LJ_2_1_6 * lj14['sigmas'][i] / 2.0
        data.nonbondedTypes.append((label, epsKcal, rMinHalf, eps14Kcal, rMin14Half))

    data.nbfixPairs = []  # [ (labelA, labelB, epsKcal, rMinFull, eps14Kcal, rMin14Full) ]
    for a in range(len(refinedTypesUsed)):
        for b in range(a + 1, len(refinedTypesUsed)):
            labelA, labelB = refinedTypesUsed[a], refinedTypesUsed[b]
            i, j = ljIndexOf[baseOfRefined[labelA]], ljIndexOf[baseOfRefined[labelB]]
            epsKcal = math.sqrt(kj_to_kcal(lj['epsilons'][i]) * kj_to_kcal(lj['epsilons'][j]))
            rMinFull = _LJ_2_1_6 * math.sqrt(lj['sigmas'][i] * lj['sigmas'][j])
            eps14Kcal = math.sqrt(kj_to_kcal(lj14['epsilons'][i]) * kj_to_kcal(lj14['epsilons'][j]))
            rMin14Full = _LJ_2_1_6 * math.sqrt(lj14['sigmas'][i] * lj14['sigmas'][j])
            data.nbfixPairs.append((labelA, labelB, epsKcal, rMinFull, eps14Kcal, rMin14Full))

    mmTerms = {type(t).__name__: t for t in mmState.mmTerms}

    def literal_orders_by_type(terms, indicesOf):
        """ (t -> representative literal refined-label tuple) is NOT
            enough once atom types are refined per-instance (see
            "_D"/"_S" above): two instances that share the same
            underlying DYFF parameter type `t` (same numeric value) can
            still end up with DIFFERENT literal refined-label tuples,
            whenever one of their atoms happens to carry a "_D" tag and
            the other doesn't -- confirmed the hard way: NAMD refused to
            load the very first generated PRM with "UNABLE TO FIND
            DIHEDRAL PARAMETERS FOR C_Tri_D C_Tri_S C_Tri_S C_Tri_S"
            because only the FIRST literal order seen for its type had
            been written. Returns EVERY distinct literal tuple actually
            observed, each mapped to one (any) `t` that produced it, so
            every literal combination NAMD could ever look up gets its
            own PRM line.
        """
        result = {}
        for term in terms:
            indices, t, isActive = indicesOf(term)
            if not isActive:
                continue
            literalOrder = tuple(data.atomTypeLabels[idx] for idx in indices)
            result.setdefault(literalOrder, t)
        return result

    # . Bonds -- exact conversion (see module docstring).
    bondState = mmTerms['HarmonicBondContainer'].__getstate__()
    data.bonds = [(i, j) for i, j, t, isActive in bondState['terms'] if isActive]
    bondParamByType = {t: (eq, kj_to_kcal(fc)) for t, (eq, fc) in enumerate(bondState['parameters'])}
    bondLiteralOrders = literal_orders_by_type(
        bondState['terms'], lambda term: ((term[0], term[1]), term[2], term[3]))
    data.bondParams = [literalOrder + bondParamByType[t] for literalOrder, t in bondLiteralOrders.items()]

    # . Angles -- curvature-matched harmonic (see module docstring).
    angleState = mmTerms['CosineAngleContainer'].__getstate__()
    angleParamByType = {}
    for t, series in enumerate(angleState['parameters']):
        angleParamByType[t] = angle_type_to_harmonic([tuple(p) for p in series])
    data.angles = [tuple(indices) for indices, t, isActive in angleState['terms'] if isActive]
    angleLiteralOrders = literal_orders_by_type(
        angleState['terms'], lambda term: (term[0], term[1], term[2]))
    data.angleParams = [literalOrder + angleParamByType[t] for literalOrder, t in angleLiteralOrders.items()]

    # . Dihedrals -- exact conversion (see module docstring).
    dihedralState = mmTerms['CosineDihedralContainer'].__getstate__()
    dihedralParamByType = {}
    for t, series in enumerate(dihedralState['parameters']):
        dihedralParamByType[t] = dihedral_type_to_charmm([tuple(p) for p in series])
    dihedralLiteralOrders = literal_orders_by_type(
        dihedralState['terms'], lambda term: (term[0], term[1], term[2]))
    data.dihedralParams = [literalOrder + dihedralParamByType[t]
                            for literalOrder, t in dihedralLiteralOrders.items()]
    data.dihedrals = [tuple(indices) for indices, t, isActive in dihedralState['terms'] if isActive]

    # . Out-of-plane / impropers -- curvature-matched harmonic, converted
    #   from DYFF's omega to CHARMM's psi via the real, per-instance
    #   geometry (see module docstring / geometric_domega_dpsi()).
    oopState = mmTerms['CosineOutOfPlaneContainer'].__getstate__()
    oopOmegaCurvatureByType = {}
    for t, series in enumerate(oopState['parameters']):
        oopOmegaCurvatureByType[t] = oop_type_curvature([tuple(p) for p in series])
    ratiosByType = {}
    data.impropers = []
    improperLiteralOrders = {}  # (labelI,labelJ,labelK,labelL) -> type index t (first one seen)
    for indices, t, isActive in oopState['terms']:
        if not isActive:
            continue
        i, j, k, l = indices
        data.impropers.append((i, j, k, l))
        ri, rj, rk, rl = data.coords[i], data.coords[j], data.coords[k], data.coords[l]
        ratio = geometric_domega_dpsi(ri, rj, rk, rl)
        if ratio is not None:
            ratiosByType.setdefault(t, []).append(ratio)
        literalOrder = (data.atomTypeLabels[i], data.atomTypeLabels[j], data.atomTypeLabels[k], data.atomTypeLabels[l])
        improperLiteralOrders.setdefault(literalOrder, (t, i, j, k, l))

    data.improperParams = []  # [ (labelI,labelJ,labelK,labelL, KpsiKcal, psi0Deg) ]
    for literalOrder, (t, i, j, k, l) in improperLiteralOrders.items():
        ratios = ratiosByType.get(t, [])
        if not ratios:
            continue
        meanRatio = sum(ratios) / len(ratios)
        _, d2eDomega2 = oopOmegaCurvatureByType[t]
        kPsiKcal = kj_to_kcal(d2eDomega2 * meanRatio * meanRatio / 2.0)
        psiRef = charmm_psi(data.coords[i], data.coords[j], data.coords[k], data.coords[l])
        psi0Deg = 0.0 if math.cos(psiRef) > 0.0 else 180.0
        data.improperParams.append(literalOrder + (kPsiKcal, psi0Deg))

    # . Box (if periodic).
    data.box = None
    if system.symmetryParameters is not None:
        sp = system.symmetryParameters
        data.box = (sp.a, sp.b, sp.c, sp.alpha, sp.beta, sp.gamma)

    data.electrostaticScale14 = system.mmModel.electrostaticScale14

    return data


#=====================================================================================
# . Writers.
#=====================================================================================
def _psf_int_lines(pairs_flat, perLine):
    """ pairs_flat: a flat list of 1-based ints, already in the right
        grouping order. Yields text lines, `perLine` GROUPS (not ints)
        per line -- e.g. perLine=4 for bonds (2 ints/group -> 8 ints
        per line), perLine=2 for dihedrals (4 ints/group -> 8 ints per
        line).
    """
    line = []
    for value in pairs_flat:
        line.append("{:8d}".format(value))
        if len(line) == perLine:
            yield "".join(line)
            line = []
    if line:
        yield "".join(line)


def write_psf(path, data):
    lines = []
    lines.append("PSF")
    lines.append("")
    lines.append("{:8d} !NTITLE".format(1))
    lines.append("* Generated by dyff_to_charmm.py -- DYFF (UFF) system exported for NAMD.")
    lines.append("")

    lines.append("{:8d} !NATOM".format(data.n))
    for i in range(data.n):
        resName, resNumber, atomName = data.residues[i]
        # . Every field explicitly space-separated (not just column-
        #   padded) -- a value that exactly fills its nominal width
        #   (e.g. the 5-character sanitized type "O_Tet") would
        #   otherwise glue onto the next field with zero space between
        #   them, corrupting NAMD's whitespace-tokenized PSF reader.
        lines.append("{:8d} {:<4s} {:<4d} {:<6s} {:<6s} {:<6s} {:12.6f} {:10.4f} {:6d}".format(
            i + 1, "SYS", resNumber, resName, atomName,
            data.atomTypeLabels[i], data.charges[i], data.masses[i], 0))
    lines.append("")

    bondInts = []
    for i, j in data.bonds:
        bondInts += [i + 1, j + 1]
    lines.append("{:8d} !NBOND: bonds".format(len(data.bonds)))
    lines.extend(_psf_int_lines(bondInts, 4))
    lines.append("")

    angleInts = []
    for i, j, k in data.angles:
        angleInts += [i + 1, j + 1, k + 1]
    lines.append("{:8d} !NTHETA: angles".format(len(data.angles)))
    lines.extend(_psf_int_lines(angleInts, 3))
    lines.append("")

    dihedralInts = []
    for i, j, k, l in data.dihedrals:
        dihedralInts += [i + 1, j + 1, k + 1, l + 1]
    lines.append("{:8d} !NPHI: dihedrals".format(len(data.dihedrals)))
    lines.extend(_psf_int_lines(dihedralInts, 2))
    lines.append("")

    improperInts = []
    for i, j, k, l in data.impropers:
        improperInts += [i + 1, j + 1, k + 1, l + 1]
    lines.append("{:8d} !NIMPHI: impropers".format(len(data.impropers)))
    lines.extend(_psf_int_lines(improperInts, 2))
    lines.append("")

    lines.append("{:8d} !NDON: donors".format(0))
    lines.append("")
    lines.append("{:8d} !NACC: acceptors".format(0))
    lines.append("")
    lines.append("{:8d} !NNB".format(0))
    lines.extend(_psf_int_lines([0] * data.n, 8))
    lines.append("")

    lines.append("{:8d}{:8d} !NGRP NST2".format(1, 0))
    lines.append("{:8d}{:8d}{:8d}".format(0, 0, 0))
    lines.append("")

    lines.append("{:8d}{:8d} !NUMLP NUMLPH".format(0, 0))

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_prm(path, data):
    lines = []
    lines.append("* Generated by dyff_to_charmm.py -- DYFF (UFF) parameters exported for NAMD.")
    lines.append("* Angle and improper terms are a curvature-matched HARMONIC approximation")
    lines.append("* to DYFF's general Fourier bends -- see dyff_to_charmm.py's module docstring.")
    lines.append("*")
    lines.append("")

    lines.append("BONDS")
    for labelI, labelJ, b0, kb in data.bondParams:
        lines.append("{:<6s} {:<6s} {:12.4f} {:12.6f}".format(labelI, labelJ, kb, b0))
    lines.append("")

    lines.append("ANGLES")
    for labelI, labelJ, labelK, theta0, ktheta in data.angleParams:
        lines.append("{:<6s} {:<6s} {:<6s} {:12.4f} {:12.4f}".format(labelI, labelJ, labelK, ktheta, theta0))
    lines.append("")

    lines.append("DIHEDRALS")
    for labelI, labelJ, labelK, labelL, kchi, n, delta in data.dihedralParams:
        lines.append("{:<6s} {:<6s} {:<6s} {:<6s} {:12.6f} {:4d} {:8.2f}".format(
            labelI, labelJ, labelK, labelL, kchi, n, delta))
    lines.append("")

    lines.append("IMPROPERS")
    for labelI, labelJ, labelK, labelL, kpsi, psi0 in data.improperParams:
        lines.append("{:<6s} {:<6s} {:<6s} {:<6s} {:12.6f} {:4d} {:8.2f}".format(
            labelI, labelJ, labelK, labelL, kpsi, 0, psi0))
    lines.append("")

    lines.append("NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfshift -")
    lines.append("cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5")
    lines.append("!atom  ignored    epsilon      Rmin/2      ignored    eps,1-4     Rmin/2,1-4")
    for label, eps, rMinHalf, eps14, rMin14Half in data.nonbondedTypes:
        lines.append("{:<6s} {:8.4f} {:12.6f} {:12.6f} {:8.4f} {:12.6f} {:12.6f}".format(
            label, 0.0, -eps, rMinHalf, 0.0, -eps14, rMin14Half))
    lines.append("")

    lines.append("NBFIX")
    lines.append("!atom1  atom2   Emin       Rmin      Emin,1-4    Rmin,1-4")
    for labelA, labelB, eps, rMin, eps14, rMin14 in data.nbfixPairs:
        lines.append("{:<6s} {:<6s} {:12.6f} {:12.6f} {:12.6f} {:12.6f}".format(
            labelA, labelB, -eps, rMin, -eps14, rMin14))
    lines.append("")
    lines.append("END")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_pdb(path, data):
    """ Strict-column PDB v3 ATOM records (name cols 13-16, altLoc 17,
        resName 18-20, chainID 22, resSeq 23-26, ... -- NAMD's own coord
        reading only needs matching ATOM COUNT/ORDER against the PSF,
        but real fixed columns keep this readable by any other tool
        too, e.g. VMD).
    """
    lines = []
    for i in range(data.n):
        resName, resNumber, atomName = data.residues[i]
        x, y, z = data.coords[i]
        lines.append(
            "ATOM  {:>5d} {:<4s}{:1s}{:>3s} {:1s}{:>4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
                (i + 1) % 100000, atomName[:4], "", resName[:3], "A", resNumber % 10000, "",
                x, y, z, 1.00, 0.00, data.elementSymbols[i], ""))
    lines.append("END")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


#=====================================================================================
# . Orchestration.
#=====================================================================================
def convert(yaml_path, output_prefix, pdynamo3_home=None):
    system = load_system(yaml_path, pdynamo3_home=pdynamo3_home)
    data = Extract(system)
    write_psf(output_prefix + ".psf", data)
    write_prm(output_prefix + ".prm", data)
    write_pdb(output_prefix + ".pdb", data)
    report = {
        "natoms": data.n,
        "nbonds": len(data.bonds),
        "nangles": len(data.angles),
        "ndihedrals": len(data.dihedrals),
        "nimpropers": len(data.impropers),
        "nAtomTypes": len(data.nonbondedTypes),
        "nNbfixPairs": len(data.nbfixPairs),
        "electrostaticScale14": data.electrostaticScale14,
        "box": data.box,
    }
    return report


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: python3 dyff_to_charmm.py <system.yaml> <output_prefix> [pdynamo3_home]\n")
        sys.exit(1)
    yaml_path = sys.argv[1]
    output_prefix = sys.argv[2]
    pdynamo3_home = sys.argv[3] if len(sys.argv) > 3 else None
    report = convert(yaml_path, output_prefix, pdynamo3_home=pdynamo3_home)
    for key, value in report.items():
        print("{}: {}".format(key, value))
