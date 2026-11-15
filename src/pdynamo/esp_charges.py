"""
esp_charges
===========

ESP and RESP atomic charges from a pDynamo HF/DFT wavefunction.

pDynamo provides the pieces that need the wavefunction: the Lebedev grids
used to build the van der Waals surfaces and the electrostatic potential of
the QC density at the grid points (QCModelBase.GridPointPotentials). The
fitting itself is done here, in numpy, because pSimulation.ESPChargeFitting
has three limits for force-field work:

  - its RESP iterator adds a*q/sqrt(q^2+b^2) to the diagonal, while the
    derivative of the hyperbolic restraint gives a/sqrt(q^2+b^2);
  - it restrains every atom (AMBER RESP leaves hydrogens free);
  - it has no equivalencing and no second stage.

The fit below follows Bayly et al., J. Phys. Chem. 97, 10269 (1993), as in
AMBER's resp/respgen:

  ESP       plain least squares, total charge constrained (Merz-Kollman
            style surfaces at 1.4, 1.6, 1.8 and 2.0 x vdW radii);
  RESP-1    a = 0.0005, b = 0.1 on heavy atoms; topologically equivalent
            atoms share one charge, except CH2/CH3 carbons and hydrogens;
  RESP-2    a = 0.001; only CH2/CH3 carbons and their hydrogens are refitted
            (equivalent hydrogens share one charge), every other atom keeps
            its RESP-1 charge.

Everything is in atomic units (bohr, hartree/e), as in the AMBER programs.

    result = compute_esp_charges(Z, xyz, bonds, charge=0, multiplicity=1,
                                 level="HF/6-31G*")
"""

import math
import time

import numpy as np

# . QC levels offered (functional, orbital basis). HF/6-31G* is the level of
#   the AMBER force fields (GAFF, ff14SB residue charges).
LEVELS = {
    "HF/6-31G*"      : ("hf"   , "6-31g_st"),
    "B3LYP/6-31G*"   : ("b3lyp", "6-31g_st"),
    "HF/def2-SVP"    : ("hf"   , "def2-svp"),
    "B3LYP/def2-SVP" : ("b3lyp", "def2-svp"),
    "PBE/def2-SVP"   : ("pbe"  , "def2-svp"),
}
# . GFN2-xTB levels: the potential of the xTB density on xtb's own vdW grid
#   (xtb --esp). Fast, but not equivalent to the HF/6-31G* RESP of the AMBER
#   force fields: differences of ~0.05 e on average and up to ~0.2 e per atom
#   (tests of 2026-10-23). GFN1 was worse (up to 0.4 e) and is not offered.
XTB_LEVELS = {"GFN2-xTB": None, "GFN2-xTB/ALPB(water)": "water"}
for _name, _solvent in XTB_LEVELS.items():
    LEVELS[_name] = ("xtb", _solvent)
DEFAULT_LEVEL = "HF/6-31G*"

SCALING_FACTORS = (1.4, 1.6, 1.8, 2.0)
LEBEDEV_ORDER   = 21            # . same as pSimulation.ESPChargeFitting

RESP_A1 = 0.0005
RESP_A2 = 0.001
RESP_B  = 0.1

BOHR_PER_ANGSTROM = 1.0 / 0.52917721092
DEBYE_PER_AU      = 2.541746473


# =============================================================================
#  Topology
# =============================================================================
_COVALENT_RADII = {1: 0.31, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57,
                   14: 1.11, 15: 1.07, 16: 1.05, 17: 1.02, 35: 1.20, 53: 1.39}


def bonds_from_distances(Z, xyz, scale=1.25):
    """Bonds (i, j) from covalent radii, for structures without connectivity."""
    xyz = np.asarray(xyz, dtype=float)
    r   = np.array([_COVALENT_RADII.get(int(z), 1.5) for z in Z])
    d   = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    cut = scale * (r[:, None] + r[None, :])
    i, j = np.where(np.triu((d < cut) & (d > 0.1), 1))
    return list(zip(i.tolist(), j.tolist()))


def _neighbours(n, bonds):
    nb = [[] for _ in range(n)]
    for i, j in bonds:
        nb[i].append(j)
        nb[j].append(i)
    return nb


def topological_classes(Z, bonds):
    """Atoms that are equivalent by topology get the same class number.

    Iterative refinement of (element, sorted neighbour classes) - the
    Morgan / Weisfeiler-Lehman scheme - until the partition stops changing.
    """
    n  = len(Z)
    nb = _neighbours(n, bonds)
    labels = [(int(Z[i]), len(nb[i])) for i in range(n)]
    classes = _renumber(labels)
    for _ in range(n + 1):
        labels = [(classes[i], tuple(sorted(classes[j] for j in nb[i]))) for i in range(n)]
        new = _renumber(labels)
        if len(set(new)) == len(set(classes)):
            classes = new
            break
        classes = new
    return classes


def _renumber(labels):
    table = {}
    for lab in sorted(set(labels)):
        table[lab] = len(table)
    return [table[lab] for lab in labels]


def methyl_methylene_groups(Z, bonds):
    """Carbons with four neighbours and at least two hydrogens: {carbon: [H, ...]}."""
    nb = _neighbours(len(Z), bonds)
    groups = {}
    for i, z in enumerate(Z):
        if int(z) != 6 or len(nb[i]) != 4:
            continue
        hs = [j for j in nb[i] if int(Z[j]) == 1]
        if len(hs) >= 2:
            groups[i] = hs
    return groups


# =============================================================================
#  Grid
# =============================================================================
def vdw_surface_points(Z, xyz, scaling_factors=SCALING_FACTORS, order=LEBEDEV_ORDER):
    """Points on scaled van der Waals surfaces (angstroms).

    The same construction as pSimulation.ESPChargeFitting.GenerateVanDerWaalsSurface
    (pDynamo radii and Lebedev grid), vectorized: a point of atom i on the
    surface at factor f is kept when it is outside f*radius of every other atom.
    """
    from pScientific import PeriodicTable
    from pScientific.Quadrature import LebedevLaikovGrid_GetGridPoints

    xyz   = np.asarray(xyz, dtype=float)
    radii = np.array([PeriodicTable[int(z)].vdwRadius for z in Z])
    basic, _ = LebedevLaikovGrid_GetGridPoints(order)
    unit  = np.array([[basic[p, 0], basic[p, 1], basic[p, 2]] for p in range(basic.rows)])
    keep  = []
    for f in scaling_factors:
        r2 = (f * radii) ** 2
        for i in range(len(Z)):
            pts = xyz[i] + f * radii[i] * unit
            d2  = ((pts[:, None, :] - xyz[None, :, :]) ** 2).sum(axis=2)
            d2[:, i] = np.inf
            ok  = (d2 > r2[None, :]).all(axis=1)
            keep.append(pts[ok])
    return np.vstack(keep) if keep else np.zeros((0, 3))


# =============================================================================
#  Fitting
# =============================================================================
def _fit(R, phi, total_charge, free_groups, fixed, restrained, a, b,
         max_iterations=500, tolerance=1.0e-8):
    """Constrained (R)ESP fit.

    R            (n, g) matrix of 1/r_ig (bohr^-1)
    phi          (g,)   reference potential (a.u.)
    free_groups  list of lists of atoms; atoms in one group share one charge
    fixed        {atom: charge} held at a value
    restrained   boolean mask of atoms that carry the hyperbolic restraint
    """
    n = R.shape[0]
    m = len(free_groups)
    T = np.zeros((n, m))
    for k, group in enumerate(free_groups):
        T[group, k] = 1.0
    c = np.zeros(n)
    for i, q in fixed.items():
        c[i] = q
    RRt  = R @ R.T
    Rphi = R @ phi
    ones = T.sum(axis=0)
    mask = np.asarray(restrained, dtype=bool).copy()
    for i in fixed:
        mask[i] = False

    q = c.copy()
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        w = np.zeros(n)
        if a > 0.0:
            w[mask] = a / np.sqrt(q[mask] ** 2 + b * b)
        H   = RRt + np.diag(w)
        K   = np.zeros((m + 1, m + 1))
        K[:m, :m] = T.T @ H @ T
        K[:m, m]  = ones
        K[m, :m]  = ones
        rhs = np.zeros(m + 1)
        rhs[:m] = T.T @ (Rphi - H @ c)
        rhs[m]  = total_charge - c.sum()
        x, *_ = np.linalg.lstsq(K, rhs, rcond=None)
        q_new = T @ x[:m] + c
        dq = np.abs(q_new - q).max()
        q = q_new
        if a == 0.0 or dq < tolerance:
            break
    residual = phi - R.T @ q
    rrms = math.sqrt(float(residual @ residual) / float(phi @ phi))
    rms  = math.sqrt(float(residual @ residual) / len(phi))
    return q, rrms, rms, iterations


def fit_charges(Z, xyz, bonds, phi, points, total_charge):
    """ESP, RESP-1 and RESP-2 charges from potentials at the grid points.

    xyz and points in angstroms, phi in a.u. Returns a dict with the three
    charge sets and their fit statistics.
    """
    Z   = [int(z) for z in Z]
    n   = len(Z)
    xyz = np.asarray(xyz, dtype=float) * BOHR_PER_ANGSTROM
    pts = np.asarray(points, dtype=float) * BOHR_PER_ANGSTROM
    R   = 1.0 / np.linalg.norm(pts[None, :, :] - xyz[:, None, :], axis=2)
    phi = np.asarray(phi, dtype=float)
    heavy = np.array([z > 1 for z in Z])

    classes = topological_classes(Z, bonds)
    groups  = methyl_methylene_groups(Z, bonds)
    in_ch   = set(groups)
    for hs in groups.values():
        in_ch.update(hs)

    out = {"classes": classes, "methyl_methylene": groups}

    # . plain ESP: every atom free, no equivalencing
    q, rrms, rms, _ = _fit(R, phi, total_charge, [[i] for i in range(n)], {}, heavy, 0.0, RESP_B)
    out["esp"] = {"charges": q, "rrms": rrms, "rms": rms}

    # . RESP-1: equivalent atoms share a charge, except CH2/CH3 carbons and hydrogens
    by_class = {}
    stage1 = []
    for i in range(n):
        if i in in_ch:
            stage1.append([i])
        else:
            by_class.setdefault(classes[i], []).append(i)
    stage1.extend(by_class.values())
    q1, rrms, rms, it1 = _fit(R, phi, total_charge, stage1, {}, heavy, RESP_A1, RESP_B)
    out["resp1"] = {"charges": q1, "rrms": rrms, "rms": rms, "iterations": it1}

    # . RESP-2: refit CH2/CH3 groups only, equivalent ones together
    if groups:
        by_class = {}
        for i in sorted(in_ch):
            by_class.setdefault(classes[i], []).append(i)
        stage2 = list(by_class.values())
        fixed  = {i: float(q1[i]) for i in range(n) if i not in in_ch}
        q2, rrms, rms, it2 = _fit(R, phi, total_charge, stage2, fixed, heavy, RESP_A2, RESP_B)
    else:
        q2, it2 = q1.copy(), 0
    out["resp2"] = {"charges": q2, "rrms": rrms, "rms": rms, "iterations": it2}
    out["points"] = len(phi)
    return out


def dipole_from_charges(xyz, charges):
    """Dipole (debye) of point charges about the centre of charge-neutral origin."""
    xyz = np.asarray(xyz, dtype=float) * BOHR_PER_ANGSTROM
    mu  = (np.asarray(charges)[:, None] * (xyz - xyz.mean(axis=0))).sum(axis=0)
    return np.linalg.norm(mu) * DEBYE_PER_AU


# =============================================================================
#  pDynamo calculation
# =============================================================================
def qc_potentials(Z, xyz, points, charge=0, multiplicity=1, level=DEFAULT_LEVEL,
                  max_scf_iterations=250, log=None):
    """SCF with pDynamo and the potential of the QC density + nuclei at points.

    Z, xyz (angstroms) of the molecule; points in angstroms.
    Returns (energy_kJ_mol, phi (a.u.), qc_dipole_debye).
    """
    import pBabel  # noqa: F401  (import order: pBabel before pMolecule)
    from pMolecule import System
    from pMolecule.QCModel import ElectronicState
    from pMolecule.QCModel import QCModelDFT, DIISSCFConverger
    from pScientific.Geometry3 import Coordinates3

    functional, basis = LEVELS[level]
    system = System.FromAtoms(list(int(z) for z in Z))
    crd = Coordinates3.WithExtent(len(Z))
    for i, (x, y, z) in enumerate(xyz):
        crd[i, 0], crd[i, 1], crd[i, 2] = float(x), float(y), float(z)
    system.coordinates3    = crd
    system.electronicState = ElectronicState.WithOptions(charge=int(charge), multiplicity=int(multiplicity))
    converger = DIISSCFConverger.WithOptions(densityTolerance=1.0e-8, maximumIterations=int(max_scf_iterations))
    # . pDynamo keeps every two-electron integral in memory (~N^4): the default 2 GB limit
    #   stops at ~22 atoms for HF/6-31G*. Allow up to half of the free memory; larger
    #   molecules need ORCA (orca_potentials), which compute_esp_charges prefers.
    options = dict(converger=converger, functional=functional, orbitalBasis=basis,
                   maximumMemory=max(2.0, 0.5 * _available_memory_gb()))
    if int(multiplicity) != 1:
        options["isSpinRestricted"] = False
    system.DefineQCModel(QCModelDFT.WithOptions(**options))
    energy = system.Energy(log=log)

    grid = Coordinates3.WithExtent(len(points))
    for p, (x, y, z) in enumerate(points):
        grid[p, 0] = x * BOHR_PER_ANGSTROM
        grid[p, 1] = y * BOHR_PER_ANGSTROM
        grid[p, 2] = z * BOHR_PER_ANGSTROM
    values = system.qcModel.GridPointPotentials(system, grid, includeNuclear=True)
    phi = np.array([values[p] for p in range(len(points))])

    try:
        mu = system.qcModel.DipoleMoment(system)
        # . pDynamo returns the dipole in debye
        qc_dipole = math.sqrt(mu[0] ** 2 + mu[1] ** 2 + mu[2] ** 2)
    except Exception:
        qc_dipole = None
    return energy, phi, qc_dipole


def find_xtb():
    """xtb executable: $PDYNAMO3_XTBCOMMAND, PATH, $XTBHOME/bin/xtb or ~/programs/xtb-*/bin/xtb."""
    import glob
    import os
    import shutil
    candidates = []
    if os.getenv("PDYNAMO3_XTBCOMMAND"):
        candidates.append(os.getenv("PDYNAMO3_XTBCOMMAND"))
    if shutil.which("xtb"):
        candidates.append(shutil.which("xtb"))
    if os.getenv("XTBHOME"):
        candidates.append(os.path.join(os.getenv("XTBHOME"), "bin", "xtb"))
    candidates += sorted(glob.glob(os.path.expanduser("~/programs/xtb-*/bin/xtb")), reverse=True)
    for path in candidates:
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return None


def xtb_potentials(Z, xyz, charge=0, multiplicity=1, solvent=None, log_dir=None):
    """GFN2-xTB single point with --esp. Returns (energy kJ/mol, points A, phi a.u., dipole D)."""
    import os
    import re
    import subprocess
    import tempfile
    from pScientific import PeriodicTable
    xtb = find_xtb()
    if xtb is None:
        raise ValueError("xtb was not found (set PDYNAMO3_XTBCOMMAND or put xtb in the PATH).")
    work = tempfile.mkdtemp(prefix="easyhybrid_xtb_esp_")
    with open(os.path.join(work, "mol.xyz"), "w") as handle:
        handle.write("{}\n\n".format(len(Z)))
        for z, (x, y, w) in zip(Z, xyz):
            handle.write("{} {:.8f} {:.8f} {:.8f}\n".format(PeriodicTable.Symbol(int(z)), x, y, w))
    command = [xtb, "mol.xyz", "--gfn", "2", "--chrg", str(int(charge)), "--uhf", str(int(multiplicity) - 1), "--esp"]
    if solvent:
        command += ["--alpb", solvent]
    env = dict(os.environ)
    env.setdefault("OMP_STACKSIZE", "1G")
    run = subprocess.run(command, cwd=work, capture_output=True, text=True, env=env)
    from pdynamo import param_logs
    param_logs.keep(work, log_dir, ["mol.xyz"], {"xtb_esp.out": "$ " + " ".join(command[1:]) + "\n" + run.stdout + run.stderr})
    path = os.path.join(work, "xtb_esp.dat")
    if run.returncode != 0 or not os.path.exists(path):
        raise ValueError("xtb failed:\n" + (run.stdout + run.stderr)[-1500:])
    data = np.loadtxt(path)
    points = data[:, :3] / BOHR_PER_ANGSTROM
    phi = data[:, 3]
    energy = float("nan")
    m = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", run.stdout)
    if m:
        energy = float(m.group(1)) * 2625.4996394799
    dipole = None
    lines = run.stdout.splitlines()
    for k, line in enumerate(lines):
        if line.strip().startswith("molecular dipole:"):
            for follow in lines[k + 1:k + 5]:
                if follow.strip().startswith("full:"):
                    dipole = float(follow.split()[-1])
                    break
    return energy, points, phi, dipole


def _available_memory_gb():
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024.0 ** 2
    except Exception:
        pass
    return 4.0


ORCA_KEYWORDS = {"HF/6-31G*": "HF 6-31G*", "B3LYP/6-31G*": "B3LYP 6-31G*", "HF/def2-SVP": "HF def2-SVP",
                 "B3LYP/def2-SVP": "B3LYP def2-SVP", "PBE/def2-SVP": "PBE def2-SVP"}


def find_orca():
    """ORCA executable with orca_vpot beside it: $PDYNAMO3_ORCACOMMAND or PATH."""
    import os
    import shutil
    for path in (os.getenv("PDYNAMO3_ORCACOMMAND"), shutil.which("orca")):
        if path and os.path.isfile(path) and os.path.isfile(os.path.join(os.path.dirname(path), "orca_vpot")):
            return path
    return None


def orca_potentials(Z, xyz, points, charge=0, multiplicity=1, level=DEFAULT_LEVEL, log_dir=None, nprocs=None):
    """SCF with ORCA and the potential at the points (orca_vpot). Same values as
    qc_potentials (pDynamo) to 1e-6 a.u. for HF/6-31G* (bALA dipeptide, 2026-10-23),
    without pDynamo's in-memory integrals. Returns (energy kJ/mol, phi a.u., dipole D)."""
    import os
    import re
    import shutil
    import subprocess
    import tempfile
    from pScientific import PeriodicTable
    orca = find_orca()
    work = tempfile.mkdtemp(prefix="easyhybrid_orca_esp_")
    try:
        nprocs = int(nprocs or os.getenv("EASYHYBRID_ORCA_NPROCS", "1"))
        keywords = ORCA_KEYWORDS[level] + (" UHF" if int(multiplicity) != 1 and level.startswith("HF") else "") + " TightSCF"
        with open(os.path.join(work, "m.inp"), "w") as handle:
            handle.write("! {}\n".format(keywords))
            if nprocs > 1:
                handle.write("%pal nprocs {} end\n".format(nprocs))
            handle.write("* xyz {} {}\n".format(int(charge), int(multiplicity)))
            for z, (x, y, w) in zip(Z, xyz):
                handle.write("{} {:.8f} {:.8f} {:.8f}\n".format(PeriodicTable.Symbol(int(z)), x, y, w))
            handle.write("*\n")
        run = subprocess.run([orca, "m.inp"], cwd=work, capture_output=True, text=True)
        from pdynamo import param_logs
        param_logs.keep(work, log_dir, ["m.inp"], {"orca.out": run.stdout + run.stderr})
        if "ORCA TERMINATED NORMALLY" not in run.stdout:
            raise ValueError("ORCA failed:\n" + (run.stdout + run.stderr)[-1500:])
        with open(os.path.join(work, "p.xyz"), "w") as handle:
            handle.write("{}\n".format(len(points)))
            for p in np.asarray(points) * BOHR_PER_ANGSTROM:
                handle.write("{:.8f} {:.8f} {:.8f}\n".format(*p))
        vpot = os.path.join(os.path.dirname(orca), "orca_vpot")
        subprocess.run([vpot, "m.gbw", "m.scfp", "p.xyz", "v.out"], cwd=work, capture_output=True, text=True)
        if not os.path.exists(os.path.join(work, "v.out")):
            raise ValueError("orca_vpot gave no potential.")
        phi = np.loadtxt(os.path.join(work, "v.out"), skiprows=1)[:, 3]
        energy = float("nan")
        m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", run.stdout)
        if m:
            energy = float(m[-1]) * 2625.4996394799
        dipole = None
        m = re.findall(r"Magnitude \(Debye\)\s*:\s*(-?\d+\.\d+)", run.stdout)
        if m:
            dipole = float(m[-1])
        return energy, phi, dipole
    finally:
        shutil.rmtree(work, ignore_errors=True)


def compute_esp_charges(Z, xyz, bonds=None, charge=0, multiplicity=1,
                        level=DEFAULT_LEVEL, log=None, engine=None, log_dir=None, nprocs=None):
    """Full calculation: grid, SCF, potentials and the three fits.

    Z, xyz (angstroms) and bonds (local indices) of one or more whole
    molecules. Returns a plain dict (it is sent back from a subprocess).
    log_dir: the original xtb/ORCA inputs and outputs are kept there;
    nprocs: ORCA cores (default $EASYHYBRID_ORCA_NPROCS or 1).
    """
    t0  = time.time()
    xyz = np.asarray(xyz, dtype=float)
    if bonds is None:
        bonds = bonds_from_distances(Z, xyz)
    if LEVELS[level][0] == "xtb":
        energy, points, phi, qc_dipole = xtb_potentials(Z, xyz, charge, multiplicity, LEVELS[level][1], log_dir=log_dir)
    else:
        points = vdw_surface_points(Z, xyz)
        engine = engine or ("orca" if find_orca() else "pdynamo")
        if engine == "orca":
            energy, phi, qc_dipole = orca_potentials(Z, xyz, points, charge, multiplicity, level, log_dir=log_dir, nprocs=nprocs)
        else:
            energy, phi, qc_dipole = qc_potentials(Z, xyz, points, charge, multiplicity, level, log=log)
    res = fit_charges(Z, xyz, bonds, phi, points, float(charge))
    for key in ("esp", "resp1", "resp2"):
        res[key]["charges"] = [float(q) for q in res[key]["charges"]]
        res[key]["dipole"]  = float(dipole_from_charges(xyz, res[key]["charges"]))
    res.update({"engine": "xtb" if LEVELS[level][0] == "xtb" else engine, "level": level, "charge": int(charge), "multiplicity": int(multiplicity),
                "energy": float(energy), "qc_dipole": qc_dipole,
                "seconds": time.time() - t0, "atoms": len(Z)})
    res["methyl_methylene"] = {int(k): [int(h) for h in v] for k, v in res["methyl_methylene"].items()}
    return res


# =============================================================================
#  Systems
# =============================================================================
def molecule_from_system(system, indices):
    """(Z, xyz, bonds) of a set of atoms of a pDynamo system, bonds in local indices.

    Raises ValueError when the atoms are bonded to atoms outside the set: the
    fit is done in the gas phase and a cut covalent bond would leave a radical.
    """
    indices = sorted(int(i) for i in indices)
    local   = {g: k for k, g in enumerate(indices)}
    Z   = [system.atoms[i].atomicNumber for i in indices]
    crd = system.coordinates3
    xyz = np.array([[crd[i, 0], crd[i, 1], crd[i, 2]] for i in indices])

    bonds, cut = [], []
    connectivity = getattr(system, "connectivity", None)
    pairs = []
    if connectivity is not None and getattr(connectivity, "bonds", None) is not None:
        for bond in connectivity.bonds:
            pairs.append((bond.node1.index, bond.node2.index))
    if pairs:
        for i, j in pairs:
            a, b = i in local, j in local
            if a and b:
                bonds.append((local[i], local[j]))
            elif a or b:
                cut.append((i, j))
    else:
        bonds = bonds_from_distances(Z, xyz)
    if cut:
        i, j = cut[0]
        raise ValueError("The selection cuts {} covalent bond(s) (for example atoms {} and {}). "
                         "Select whole molecules (selecting by molecule) - the fit is done in the "
                         "gas phase.".format(len(cut), i, j))
    return Z, xyz, bonds


def suggested_charge(system, indices):
    """Net charge from the MM charges of the atoms, rounded; 0 without an MM model."""
    try:
        charges = system.mmState.charges
        return int(round(sum(charges[i] for i in indices)))
    except Exception:
        return 0


def apply_to_mm_model(system, indices, charges):
    """Writes charges into system.mmState.charges. Returns (old_sum, new_sum)."""
    mm = getattr(system, "mmState", None)
    if mm is None:
        raise ValueError("The system has no MM model.")
    old = sum(mm.charges[i] for i in indices)
    for i, q in zip(sorted(indices), charges):
        mm.charges[i] = float(q)
    return old, sum(charges)


# =============================================================================
#  Command line (the GUI runs the calculation in a separate Python process)
# =============================================================================
def run_job_file(job_path, result_path):
    """Reads a JSON job {Z, xyz, bonds, charge, multiplicity, level} and writes
    the JSON result, or {"error": text}."""
    import json
    import traceback
    with open(job_path) as handle:
        job = json.load(handle)
    try:
        result = compute_esp_charges(job["Z"], job["xyz"], [tuple(b) for b in job["bonds"]],
                                     job["charge"], job["multiplicity"], job["level"])
        result["methyl_methylene"] = {str(k): v for k, v in result["methyl_methylene"].items()}
    except Exception as error:
        result = {"error": "{}\n\n{}".format(error, traceback.format_exc(limit=3))}
    with open(result_path + ".tmp", "w") as handle:
        json.dump(result, handle)
    import os
    os.replace(result_path + ".tmp", result_path)


if __name__ == "__main__":
    import sys
    run_job_file(sys.argv[1], sys.argv[2])
