"""
orbital_analysis
================

Molecular orbitals of a pDynamo QC system and their composition in atomic
orbitals (the LCAO), for the Orbitals and Surfaces window. No GTK.

    wf = compute_wavefunction(system)          # runs system.Energy()
    rows = orbital_composition(wf, k)          # [(atom, ao label, c, population), ...]
    text = orbital_character(wf, k)            # "O3 2p 71%, C2 2p 12%"

The coefficients C (basis functions x orbitals) come from scratch.orbitalsP,
already in the non-orthogonal atomic basis (QCModelBase.GetNonOrthogonalOrbitals).

Population of basis function mu in orbital k:
  - DFT/HF (non-orthogonal Gaussian basis, scratch.overlapMatrix):
    Mulliken, P = C_mu,k (S C)_mu,k; sums to 1 over mu.
  - NDDO methods (MNDO, AM1, PM3, PM6, RM1...): the basis is orthogonal by
    assumption (ZDO), so P = C_mu,k^2; sums to 1.
Mulliken populations of single functions can be slightly negative or above 1.

pDynamo's labels are real spherical harmonics: p0 = pz, p1c = px, p1s = py,
d0 = dz2, d1c = dxz, d1s = dyz, d2c = dx2-y2, d2s = dxy.
"""

import numpy as np

HARTREE_TO_EV = 27.211386245988

_AO_NAMES = {
    "p0": "pz", "p1c": "px", "p1s": "py",
    "d0": "dz2", "d1c": "dxz", "d1s": "dyz", "d2c": "dx2-y2", "d2s": "dxy",
    "f0": "fz3", "f1c": "fxz2", "f1s": "fyz2", "f2c": "fz(x2-y2)", "f2s": "fxyz", "f3c": "fx(x2-3y2)",
    "f3s": "fy(3x2-y2)",
}


def ao_name(label):
    """ '2p1c' -> '2px', '1s' -> '1s', '4d2s' -> '4dxy'. """
    label = str(label)
    k = 0
    while k < len(label) and label[k].isdigit():
        k += 1
    n, rest = label[:k], label[k:]
    return n + _AO_NAMES.get(rest, rest)


def ao_shell(label):
    """ '2p1c' -> '2p' (shell: principal number and angular momentum). """
    label = str(label)
    k = 0
    while k < len(label) and label[k].isdigit():
        k += 1
    return label[:k + 1]


def _array2d(a):
    if hasattr(a, "shape") and not hasattr(a, "rows"):
        return np.asarray(a, dtype=float)
    rows, cols = a.rows, a.columns
    return np.array([[a[i, j] for j in range(cols)] for i in range(rows)], dtype=float)


def has_orbitals(system):
    """ True for pDynamo's own SCF models (QCModelBase: MNDO family, DFT, HF). """
    try:
        from pMolecule.QCModel.QCModelBase import QCModelBase
        return isinstance(system.qcModel, QCModelBase)
    except Exception:
        return False


def method_name(system):
    model = system.qcModel
    for attr in ("hamiltonian", "functional"):
        value = getattr(model, attr, None)
        if value:
            basis = getattr(model, "orbitalBasis", None)
            return "{}{}".format(str(value).upper(), "/" + str(basis) if basis and attr == "functional" else "")
    return type(model).__name__


def compute_wavefunction(system):
    """ SCF at the current coordinates and everything the window shows. """
    if not has_orbitals(system):
        raise ValueError("This QC model gives no molecular orbitals to EasyHybrid (only pDynamo's own "
                         "methods do: MNDO, AM1, PM3, PM6, RM1..., DFT and HF). Use a cube file instead.")
    system.Energy(log=None)
    scratch = system.scratch
    orbitals = scratch.orbitalsP
    C = _array2d(orbitals.orbitals)
    energies = np.array(list(orbitals.energies), dtype=float)
    occupancies = np.array(list(orbitals.occupancies), dtype=float)
    n_occupied = int(orbitals.occupancyHandler.numberOccupied)
    bases = system.qcState.orbitalBases
    centers = [int(c) for c in bases.functionCenters]
    labels = [str(l) for l in bases.functionLabels]
    S = getattr(scratch, "overlapMatrix", None)
    if S is not None:
        n = S.rows
        S = np.array([[S[i, j] for j in range(n)] for i in range(n)], dtype=float)
        if S.shape[0] != C.shape[0]:
            S = None
    if S is not None:
        populations = C * (S @ C)
        population_kind = "Mulliken"
    else:
        populations = C ** 2
        population_kind = "c^2 (orthogonal NDDO basis)"
    # . 2026-10-26: the QC centers are the pure QC atoms, in order, and then the boundary
    #   atoms, which carry a link hydrogen (pDynamo's QCModel: atomicNumbers = pure + 1 per
    #   boundary atom). qcState.qcAtoms is the sorted union: with a QC/MM boundary every label
    #   after the boundary atom was shifted by one atom.
    state = system.qcState
    pure = [int(i) for i in state.pureQCAtoms]
    boundary = [int(i) for i in (getattr(state, "boundaryAtoms", None) or [])]
    qc_atoms = pure + boundary
    numbers = [int(system.atoms[i].atomicNumber) for i in pure] + [1] * len(boundary)
    if S is None:
        # . NDDO: a valence-only basis that pDynamo labels from 1 (carbon '1s', '2p'); use the
        #   principal number of the valence shell instead (C 2s 2p, S 3s 3p 3d, Fe 4s 4p 3d).
        labels = [_valence_label(l, numbers[c]) for l, c in zip(labels, centers)]
    atom_labels = []
    for c, i in enumerate(qc_atoms):
        if c < len(pure):
            atom_labels.append("{}{}".format(_symbol(system.atoms[i].atomicNumber), i + 1))
        else:
            # . the link hydrogen that replaces the bond to boundary atom i
            atom_labels.append("HL{}".format(i + 1))
    return {
        "energies": energies,                   # hartree
        "occupancies": occupancies,
        "n_occupied": n_occupied,               # HOMO = n_occupied - 1
        "C": C,
        "populations": populations,
        "population_kind": population_kind,
        "centers": centers,                     # basis function -> QC center
        "labels": labels,                       # pDynamo function labels
        "qc_atoms": qc_atoms,                   # QC center -> system atom index
        "atom_labels": atom_labels,             # QC center -> "O3"
        "method": method_name(system),
        "n_basis": C.shape[0],
    }


def _period(z):
    for n, last in enumerate((2, 10, 18, 36, 54, 86, 118), start=1):
        if z <= last:
            return n
    return 7


def _valence_label(label, z):
    shell = ao_shell(label)
    kind = shell[-1:]
    rest = label[len(shell):]
    n = _period(z)
    transition = 21 <= z <= 30 or 39 <= z <= 48 or 57 <= z <= 80 or 89 <= z <= 112
    if kind == "d" and transition:
        n -= 1
    elif kind == "f":
        n -= 2
    return "{}{}{}".format(n, kind, rest)


def _symbol(z):
    try:
        from pScientific import PeriodicTable
        return PeriodicTable.Symbol(int(z))
    except Exception:
        return "X"


def orbital_label(wf, k):
    homo = wf["n_occupied"] - 1
    if k <= homo:
        d = k - homo
        return "HOMO" if d == 0 else "HOMO{:+d}".format(d)
    d = k - homo - 1
    return "LUMO" if d == 0 else "LUMO{:+d}".format(d)


def orbital_composition(wf, k):
    """ [(center, ao name, coefficient, population), ...], largest |population| first. """
    C, P = wf["C"], wf["populations"]
    rows = [(wf["centers"][mu], ao_name(wf["labels"][mu]), float(C[mu, k]), float(P[mu, k]))
            for mu in range(C.shape[0])]
    rows.sort(key=lambda r: -abs(r[3]))
    return rows


def atom_contributions(wf, k):
    """ {center: population} of orbital k. """
    out = {}
    P = wf["populations"]
    for mu, c in enumerate(wf["centers"]):
        out[c] = out.get(c, 0.0) + float(P[mu, k])
    return out


def shell_contributions(wf, k):
    """ {(center, shell): population}, e.g. (2, '2p'): 0.71. For Gaussian bases (DFT/HF) the
        shells of one angular momentum are summed ('p'): pDynamo numbers them in sequence
        (6-31G*: 2p inner and 3p outer valence), not by principal quantum number. """
    out = {}
    P = wf["populations"]
    gaussian = wf.get("population_kind") == "Mulliken"
    for mu, c in enumerate(wf["centers"]):
        shell = ao_shell(wf["labels"][mu])
        key = (c, shell[-1:] if gaussian else shell)
        out[key] = out.get(key, 0.0) + float(P[mu, k])
    return out


def orbital_character(wf, k, n=3, minimum=0.05):
    """ Short description: the largest atom-shell contributions, 'O3 2p 71%, C2 2p 12%'. """
    items = sorted(shell_contributions(wf, k).items(), key=lambda kv: -kv[1])
    parts = ["{} {} {:.0f}%".format(wf["atom_labels"][c], shell, 100 * p) for (c, shell), p in items[:n]
             if p >= minimum]
    return ", ".join(parts) if parts else "delocalized"
