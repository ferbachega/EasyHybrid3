"""
traj_analysis
=============

Trajectory analyses for Analysis > Trajectory Analysis (GTK-free, NumPy only).

Coordinates are in angstroms, frames arrays are (n_frames, n_atoms, 3), as
vismol_object.frames. Boxes are orthorhombic (a, b, c), as everywhere else in
EasyHybrid (vismol_object.cell_coordinates[f][-1]).

    radius_of_gyration            1  compactness per frame
    hydrogen_bonds                2  donor-H...acceptor occupancies (GROMACS criterion)
    residue_contact_map           3  contact frequencies between residues
    ca_distance_map               3  mean CA-CA distances
    pca / dccm                    4, 5  essential dynamics and cross-correlations
    rmsd_matrix / gromos_clusters 6  conformational clustering (Daura et al. 1999)
    dssp                          7  secondary structure per frame (Kabsch & Sander 1983)
    sasa                          8  solvent accessible surface (Shrake & Rupley 1973)
    mean_square_displacement      9  MSD and self-diffusion (Einstein relation)
    interaction_energy           10  Coulomb + Lennard-Jones between two groups
    side_chain_rotamers          11  chi1/chi2 populations and transitions
    block_average, autocorrelation 12  statistics of any series
    fluctuations                 B-factors, anisotropic U, average structure
    quasi_harmonic               quasi-harmonic frequencies and entropy
"""

import math

import numpy as np

# =============================================================================
#  Constants and element data
# =============================================================================
COULOMB_KJ   = 1389.35457      # . kJ mol^-1 A e^-2  (e^2 / 4 pi eps0)
KB_KJ        = 0.0083144626    # . kJ mol^-1 K^-1
HBAR_JS      = 1.054571817e-34
KB_JK        = 1.380649e-23
AMU_KG       = 1.66053906660e-27
AVOGADRO     = 6.02214076e23
C_CM_S       = 2.99792458e10

MASSES = {"H": 1.008, "D": 2.014, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "NA": 22.990,
          "MG": 24.305, "P": 30.974, "S": 32.06, "CL": 35.45, "K": 39.098, "CA": 40.078, "MN": 54.938,
          "FE": 55.845, "CO": 58.933, "NI": 58.693, "CU": 63.546, "ZN": 65.38, "BR": 79.904,
          "SE": 78.971, "I": 126.904, "CD": 112.414, "HG": 200.592, "LI": 6.94, "B": 10.81, "SI": 28.085}

# . Bondi radii (A) for SASA; Shrake-Rupley with a 1.4 A probe
VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80, "CL": 1.75,
             "BR": 1.85, "I": 1.98, "SE": 1.90, "NA": 2.27, "K": 2.75, "MG": 1.73, "CA": 2.31, "ZN": 1.39,
             "FE": 1.56, "CU": 1.40, "MN": 1.61, "LI": 1.82}


def element_of(symbol):
    s = str(symbol).strip().upper()
    return s if s in MASSES else (s[:1] if s[:1] in MASSES else s)


def masses_of(symbols):
    return np.array([MASSES.get(element_of(s), 12.0) for s in symbols])


# =============================================================================
#  Superposition
# =============================================================================
def _kabsch_rotation(P, Q, w=None):
    """Rotation R (3x3) that best maps centred P onto centred Q (P @ R ~ Q)."""
    if w is None:
        H = P.T @ Q
    else:
        H = (P * w[:, None]).T @ Q
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, 1.0, d])
    return U @ D @ Vt


def superpose(frames, reference, subset=None, weights=None):
    """Every frame superposed on reference by the atoms of subset (all atoms move)."""
    frames = np.asarray(frames, dtype=float)
    idx = np.arange(frames.shape[1]) if subset is None else np.asarray(subset)
    w = None if weights is None else np.asarray(weights, dtype=float)
    ref = np.asarray(reference, dtype=float)[idx]
    ref_c = ref.mean(axis=0) if w is None else (ref * w[:, None]).sum(0) / w.sum()
    out = np.empty_like(frames)
    for f, X in enumerate(frames):
        P = X[idx]
        pc = P.mean(axis=0) if w is None else (P * w[:, None]).sum(0) / w.sum()
        R = _kabsch_rotation(P - pc, ref - ref_c, w)
        out[f] = (X - pc) @ R + ref_c
    return out


def iterative_average(frames, subset=None, iterations=4, weights=None):
    """Superposes the frames on their own average, iterating (the average is
    recomputed after each pass). Returns (aligned_frames, average)."""
    frames = np.asarray(frames, dtype=float)
    aligned = superpose(frames, frames[0], subset, weights)
    for _ in range(iterations):
        average = aligned.mean(axis=0)
        aligned = superpose(aligned, average, subset, weights)
    return aligned, aligned.mean(axis=0)


def _batched_rmsd(A, B):
    """Minimum RMSD after optimal superposition between A (n,3) and each B[k] (m,n,3)."""
    A = A - A.mean(axis=0)
    Bc = B - B.mean(axis=1, keepdims=True)
    H = np.einsum("ni,knj->kij", A, Bc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    S[:, 2] *= d
    n = A.shape[0]
    e0 = (A * A).sum() + (Bc * Bc).sum(axis=(1, 2))
    msd = np.maximum((e0 - 2.0 * S.sum(axis=1)) / n, 0.0)
    return np.sqrt(msd)


def rmsd_matrix(frames, subset=None, progress=None):
    """Pairwise RMSD (A) after optimal superposition, (F, F)."""
    X = np.asarray(frames, dtype=float)
    if subset is not None:
        X = X[:, subset, :]
    F = X.shape[0]
    D = np.zeros((F, F))
    for i in range(F - 1):
        D[i, i + 1:] = _batched_rmsd(X[i], X[i + 1:])
        D[i + 1:, i] = D[i, i + 1:]
        if progress is not None:
            progress((i + 1) / (F - 1))
    return D


# =============================================================================
#  Neighbour search (cell lists, vectorised)
# =============================================================================
def _ragged_arange(starts, counts):
    total = int(counts.sum())
    if total == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    owner = np.repeat(np.arange(len(counts)), counts)
    offsets = np.cumsum(counts) - counts
    within = np.arange(total) - np.repeat(offsets, counts)
    return owner, np.repeat(starts, counts) + within


def pairs_within(a, b, cutoff, box=None):
    """All pairs (i, j) with |a_i - b_j| <= cutoff. Minimum image when box is given.

    Returns (i, j, distance). Uses cell lists, so it scales to solvated systems.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
    if box is not None:
        box = np.asarray(box, dtype=float)
        ncell = np.floor(box / cutoff).astype(int)
        if (ncell < 3).any():
            box_ok = False
        else:
            box_ok = True
    if box is None or not box_ok:
        if len(a) * len(b) <= 4_000_000:
            d = a[:, None, :] - b[None, :, :]
            if box is not None:
                d -= box * np.round(d / box)
            r = np.sqrt((d * d).sum(axis=2))
            i, j = np.nonzero(r <= cutoff)
            return i, j, r[i, j]
        if box is None:
            lo = np.minimum(a.min(axis=0), b.min(axis=0)) - cutoff
            hi = np.maximum(a.max(axis=0), b.max(axis=0)) + cutoff
            ncell = np.maximum(np.floor((hi - lo) / cutoff).astype(int), 1)
            cell = (hi - lo) / ncell
            periodic = False
        else:
            # . small periodic box: chunked brute force
            out_i, out_j, out_r = [], [], []
            for s in range(0, len(a), 256):
                d = a[s:s + 256, None, :] - b[None, :, :]
                d -= box * np.round(d / box)
                r = np.sqrt((d * d).sum(axis=2))
                i, j = np.nonzero(r <= cutoff)
                out_i.append(i + s); out_j.append(j); out_r.append(r[i, j])
            return np.concatenate(out_i), np.concatenate(out_j), np.concatenate(out_r)
    else:
        lo = np.zeros(3)
        cell = box / ncell
        periodic = True

    def cell_of(x):
        c = np.floor((x - lo) / cell).astype(int)
        if periodic:
            return np.mod(c, ncell)
        return np.clip(c, 0, ncell - 1)

    ca, cb = cell_of(a), cell_of(b)
    key_b = (cb[:, 0] * ncell[1] + cb[:, 1]) * ncell[2] + cb[:, 2]
    order = np.argsort(key_b, kind="stable")
    sorted_keys = key_b[order]
    n_keys = int(ncell.prod())
    starts_all = np.searchsorted(sorted_keys, np.arange(n_keys), side="left")
    ends_all = np.searchsorted(sorted_keys, np.arange(n_keys), side="right")

    out_i, out_j, out_r = [], [], []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                nc = ca + np.array([dx, dy, dz])
                if periodic:
                    nc = np.mod(nc, ncell)
                    valid = np.ones(len(a), bool)
                else:
                    valid = ((nc >= 0) & (nc < ncell)).all(axis=1)
                ia = np.nonzero(valid)[0]
                if len(ia) == 0:
                    continue
                key = (nc[ia, 0] * ncell[1] + nc[ia, 1]) * ncell[2] + nc[ia, 2]
                st, en = starts_all[key], ends_all[key]
                owner, pos = _ragged_arange(st, en - st)
                if len(pos) == 0:
                    continue
                i = ia[owner]
                j = order[pos]
                d = a[i] - b[j]
                if periodic:
                    d -= box * np.round(d / box)
                r = np.sqrt((d * d).sum(axis=1))
                keep = r <= cutoff
                out_i.append(i[keep]); out_j.append(j[keep]); out_r.append(r[keep])
    if not out_i:
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
    i = np.concatenate(out_i); j = np.concatenate(out_j); r = np.concatenate(out_r)
    if periodic and (ncell == 3).any():
        # . with 3 cells along an axis the -1 and +1 neighbours are the same cell
        pair = np.unique(np.stack([i, j], axis=1), axis=0, return_index=True)[1]
        i, j, r = i[pair], j[pair], r[pair]
    return i, j, r


# =============================================================================
#  1. Radius of gyration
# =============================================================================
def radius_of_gyration(frames, masses=None):
    X = np.asarray(frames, dtype=float)
    m = np.ones(X.shape[1]) if masses is None else np.asarray(masses, dtype=float)
    com = (X * m[None, :, None]).sum(axis=1) / m.sum()
    d2 = ((X - com[:, None, :]) ** 2).sum(axis=2)
    return np.sqrt((d2 * m[None, :]).sum(axis=1) / m.sum())


# =============================================================================
#  2. Hydrogen bonds
# =============================================================================
def neighbours_from_bonds(n, bonds):
    nb = [[] for _ in range(n)]
    for i, j in bonds:
        nb[i].append(j)
        nb[j].append(i)
    return nb


def donors_and_acceptors(elements, bonds, atoms=None):
    """Donor pairs (D, H) and acceptors among `atoms` (all when None).

    Donors: N, O, S bonded to H. Acceptors: O, S, F, and N with at most two
    neighbours (aromatic N without H, nitrile) - amide, amine and ammonium N
    are not acceptors.
    """
    elements = [element_of(e) for e in elements]
    nb = neighbours_from_bonds(len(elements), bonds)
    pool = range(len(elements)) if atoms is None else atoms
    donors, acceptors = [], []
    for i in pool:
        e = elements[i]
        if e in ("N", "O", "S"):
            for h in nb[i]:
                if elements[h] == "H":
                    donors.append((i, h))
        if e in ("O", "S", "F") or (e == "N" and len(nb[i]) <= 2):
            acceptors.append(i)
    return donors, acceptors


def hydrogen_bonds(frames, donors, acceptors, boxes=None, d_cutoff=3.5, angle_cutoff=30.0, progress=None):
    """Hydrogen bonds per frame with the GROMACS criterion:
    D...A <= d_cutoff and the H-D...A angle <= angle_cutoff (degrees).

    Returns (counts per frame, {(D, H, A): [frames]}, {(D, H, A): (mean D...A, mean angle)}).
    """
    donors = np.asarray(donors, dtype=int).reshape(-1, 2)
    acceptors = np.asarray(acceptors, dtype=int)
    cos_cut = math.cos(math.radians(angle_cutoff))
    counts, present, geo = [], {}, {}
    for f, X in enumerate(frames):
        box = None if boxes is None else boxes[f]
        D = X[donors[:, 0]]
        A = X[acceptors]
        i, j, r = pairs_within(D, A, d_cutoff, box)
        keep = donors[i, 0] != acceptors[j]
        i, j, r = i[keep], j[keep], r[keep]
        dh = X[donors[i, 1]] - D[i]
        da = A[j] - D[i]
        if box is not None:
            dh -= box * np.round(dh / box)
            da -= box * np.round(da / box)
        cosang = (dh * da).sum(axis=1) / (np.linalg.norm(dh, axis=1) * np.maximum(r, 1e-9))
        ok = cosang >= cos_cut
        n = 0
        for k in np.nonzero(ok)[0]:
            key = (int(donors[i[k], 0]), int(donors[i[k], 1]), int(acceptors[j[k]]))
            present.setdefault(key, []).append(f)
            s = geo.setdefault(key, [0.0, 0.0])
            s[0] += r[k]
            s[1] += math.degrees(math.acos(min(1.0, max(-1.0, cosang[k]))))
            n += 1
        counts.append(n)
        if progress is not None:
            progress((f + 1) / len(frames))
    geometry = {k: (v[0] / len(present[k]), v[1] / len(present[k])) for k, v in geo.items()}
    return np.array(counts), present, geometry


def lifetimes(frame_list):
    """Lengths of the runs of consecutive frames."""
    if not frame_list:
        return []
    runs, start, prev = [], frame_list[0], frame_list[0]
    for f in frame_list[1:]:
        if f != prev + 1:
            runs.append(prev - start + 1)
            start = f
        prev = f
    runs.append(prev - start + 1)
    return runs


# =============================================================================
#  3. Contacts and distance maps
# =============================================================================
def residue_contact_map(frames, residue_of_atom, n_residues, cutoff=4.5, boxes=None, progress=None):
    """Fraction of frames in which two residues have any pair of atoms closer
    than cutoff (residue_of_atom[k] = residue number 0..R-1 of atom k)."""
    res = np.asarray(residue_of_atom, dtype=int)
    freq = np.zeros((n_residues, n_residues))
    for f, X in enumerate(frames):
        box = None if boxes is None else boxes[f]
        i, j, _r = pairs_within(X, X, cutoff, box)
        ri, rj = res[i], res[j]
        keep = ri != rj
        m = np.zeros((n_residues, n_residues), dtype=bool)
        m[ri[keep], rj[keep]] = True
        freq += m
        if progress is not None:
            progress((f + 1) / len(frames))
    return freq / max(len(frames), 1)


def ca_distance_map(frames):
    """Mean and standard deviation of the distances between the atoms (e.g. CA)."""
    X = np.asarray(frames, dtype=float)
    d = np.sqrt(((X[:, :, None, :] - X[:, None, :, :]) ** 2).sum(axis=3))
    return d.mean(axis=0), d.std(axis=0)


# =============================================================================
#  4. PCA   5. DCCM
# =============================================================================
def pca(aligned_frames, masses=None):
    """Principal components of the (optionally mass-weighted) covariance.

    aligned_frames: (F, n, 3), already superposed. Returns a dict with
    eigenvalues (A^2 or amu A^2, descending), eigenvectors (3n, k),
    projections (F, k) and the variance fraction of each component.
    """
    X = np.asarray(aligned_frames, dtype=float)
    F, n, _ = X.shape
    mean = X.mean(axis=0)
    D = (X - mean).reshape(F, 3 * n)
    if masses is not None:
        D = D * np.repeat(np.sqrt(np.asarray(masses, dtype=float)), 3)[None, :]
    if F < 3 * n:
        # . the covariance has rank < F: diagonalize the F x F Gram matrix (same nonzero spectrum)
        G = D @ D.T / F
        w, U = np.linalg.eigh(G)
        order = np.argsort(w)[::-1]
        w, U = w[order], U[:, order]
        keep = w > 1e-10 * max(w[0], 1e-30)
        w, U = w[keep], U[:, keep]
        V = D.T @ U / np.sqrt(w * F)
    else:
        C = D.T @ D / F
        w, V = np.linalg.eigh(C)
        order = np.argsort(w)[::-1]
        w, V = w[order], V[:, order]
    proj = D @ V
    return {"eigenvalues": w, "eigenvectors": V, "projections": proj, "mean": mean,
            "fraction": w / w.sum() if w.sum() > 0 else w}


def mode_frames(mean, eigenvector, eigenvalue, n_frames=21, sigmas=2.0, masses=None):
    """Structures along one principal component, from -sigmas to +sigmas std."""
    v = np.asarray(eigenvector, dtype=float).reshape(-1, 3)
    if masses is not None:
        v = v / np.sqrt(np.asarray(masses, dtype=float))[:, None]
    amplitude = sigmas * math.sqrt(max(eigenvalue, 0.0))
    t = np.concatenate([np.linspace(0, 1, n_frames // 2 + 1), np.linspace(1, -1, n_frames)[1:],
                        np.linspace(-1, 0, n_frames // 2 + 1)[1:]])
    return np.array([mean + amplitude * s * v for s in t])


def dccm(aligned_frames):
    """Dynamic cross-correlation matrix C_ij = <dr_i.dr_j> / sqrt(<dr_i^2><dr_j^2>)."""
    X = np.asarray(aligned_frames, dtype=float)
    D = X - X.mean(axis=0)
    cov = np.einsum("fik,fjk->ij", D, D) / X.shape[0]
    d = np.sqrt(np.diag(cov))
    return cov / np.outer(d, d)


# =============================================================================
#  6. Clustering
# =============================================================================
def gromos_clusters(D, cutoff):
    """GROMOS clustering (Daura et al., Angew. Chem. 1999): the frame with the
    most neighbours within cutoff and its neighbours form a cluster; repeat on
    the remaining frames. Returns (labels, centres) with clusters sorted by size."""
    F = D.shape[0]
    remaining = np.ones(F, dtype=bool)
    labels = np.full(F, -1)
    centres = []
    adjacency = D <= cutoff
    while remaining.any():
        counts = (adjacency & remaining[None, :]).sum(axis=1)
        counts[~remaining] = -1
        c = int(np.argmax(counts))
        members = np.nonzero(adjacency[c] & remaining)[0]
        labels[members] = len(centres)
        centres.append(c)
        remaining[members] = False
    return labels, centres


# =============================================================================
#  7. DSSP
# =============================================================================
def _place_amide_h(N, C_prev, O_prev):
    """DSSP rule: H on N along the C=O bond of the previous residue, 1.0 A."""
    v = C_prev - O_prev
    v /= np.linalg.norm(v, axis=-1, keepdims=True)
    return N + v


def dssp(frames, backbone, chain_break_after=None, progress=None):
    """Secondary structure codes per frame (Kabsch & Sander 1983).

    backbone: list over residues of dicts with atom indices N, CA, C, O and
    optionally H. chain_break_after: residues k after which the chain is broken.
    Returns an (F, R) array of one-letter codes: H G I E B T S and '-' (coil).
    The DSSP 2.x pi-helix priority (I over H) is used.
    """
    R = len(backbone)
    breaks = set(chain_break_after or [])
    idxN = np.array([b["N"] for b in backbone]); idxCA = np.array([b["CA"] for b in backbone])
    idxC = np.array([b["C"] for b in backbone]); idxO = np.array([b["O"] for b in backbone])
    has_prev = np.array([k > 0 and (k - 1) not in breaks for k in range(R)])
    out = np.full((len(frames), R), "-", dtype="<U1")
    for f, X in enumerate(frames):
        N, CA, C, O = X[idxN], X[idxCA], X[idxC], X[idxO]
        H = np.full_like(N, np.nan)
        for k, b in enumerate(backbone):
            if b.get("H") is not None:
                H[k] = X[b["H"]]
            elif has_prev[k] and backbone[k].get("res") != "PRO":
                H[k] = _place_amide_h(N[k], C[k - 1], O[k - 1])
        # . H-bond energy E(i -> j): C=O of i accepts from N-H of j (kcal/mol)
        valid_h = ~np.isnan(H[:, 0])
        hb = np.zeros((R, R), dtype=bool)
        ii, jj, _ = pairs_within(CA, CA, 9.0)
        sel = (ii != jj) & valid_h[jj] & (np.abs(ii - jj) >= 2)
        ii, jj = ii[sel], jj[sel]
        rON = np.linalg.norm(O[ii] - N[jj], axis=1)
        rCH = np.linalg.norm(C[ii] - H[jj], axis=1)
        rOH = np.linalg.norm(O[ii] - H[jj], axis=1)
        rCN = np.linalg.norm(C[ii] - N[jj], axis=1)
        E = 0.084 * 332.0 * (1.0 / rON + 1.0 / rCH - 1.0 / rOH - 1.0 / rCN)
        E = np.maximum(E, -9.9)
        ok = E < -0.5
        hb[ii[ok], jj[ok]] = True

        def connected(a, b):
            lo, hi = min(a, b), max(a, b)
            return not any(k in breaks for k in range(lo, hi))

        # . n-turns at i: hbond(i, i+n)
        turn = {n: np.zeros(R, dtype=bool) for n in (3, 4, 5)}
        for n in (3, 4, 5):
            for i in range(R - n):
                if hb[i, i + n] and connected(i, i + n):
                    turn[n][i] = True
        ss = np.full(R, "-", dtype="<U1")
        # . bridges (i < j), with their type: parallel "p" or antiparallel "a"
        bridges = []
        for i in range(1, R - 1):
            for j in range(i + 3, R - 1):
                if (hb[i - 1, j] and hb[j, i + 1]) or (hb[j - 1, i] and hb[i, j + 1]):
                    bridges.append((i, j, "p"))
                if (hb[i, j] and hb[j, i]) or (hb[i - 1, j + 1] and hb[j - 1, i + 1]):
                    bridges.append((i, j, "a"))
        # . ladders: runs of consecutive bridges of one type
        ladders = []
        for (i, j, t) in sorted(bridges):
            for lad in ladders:
                li, lj = lad["pairs"][-1]
                if lad["type"] == t and i == li + 1 and j == (lj + 1 if t == "p" else lj - 1):
                    lad["pairs"].append((i, j))
                    break
            else:
                ladders.append({"type": t, "pairs": [(i, j)]})
        # . beta-bulges: ladders of one type joined across a gap of at most
        #   1 residue on one strand and at most 4 on the other
        strand = np.zeros(R, dtype=bool)
        bridged = np.zeros(R, dtype=bool)
        linked = [False] * len(ladders)
        for a in range(len(ladders)):
            for b in range(len(ladders)):
                A, B = ladders[a], ladders[b]
                if a == b or A["type"] != B["type"]:
                    continue
                ai_end, aj_end = A["pairs"][-1]
                bi_start, bj_start = B["pairs"][0]
                gap1 = bi_start - ai_end - 1
                gap2 = (bj_start - aj_end - 1) if A["type"] == "p" else (aj_end - bj_start - 1)
                if gap1 >= 0 and gap2 >= 0 and ((gap1 <= 1 and gap2 <= 4) or (gap1 <= 4 and gap2 <= 1)):
                    linked[a] = linked[b] = True
                    strand[ai_end:bi_start + 1] = True
                    lo, hi = sorted((aj_end, bj_start))
                    strand[lo:hi + 1] = True
        for k, lad in enumerate(ladders):
            for (i, j) in lad["pairs"]:
                bridged[i] = bridged[j] = True
                if len(lad["pairs"]) > 1 or linked[k]:
                    strand[i] = strand[j] = True
        in_ladder = strand
        # . helices (priority H > B > E > G > I > T > S, with DSSP 2 pi-helix over H)
        helix = {n: np.zeros(R, dtype=bool) for n in (3, 4, 5)}
        for n in (3, 4, 5):
            for i in range(1, R - n):
                if turn[n][i] and turn[n][i - 1]:
                    helix[n][i:i + n] = True
        for i in range(R):
            if helix[5][i]:
                ss[i] = "I"
            elif helix[4][i]:
                ss[i] = "H"
            elif in_ladder[i]:
                ss[i] = "E"
            elif bridged[i]:
                ss[i] = "B"
            elif helix[3][i]:
                ss[i] = "G"
        # . turns and bends
        for n in (3, 4, 5):
            for i in range(R - n):
                if turn[n][i]:
                    for k in range(i + 1, i + n):
                        if ss[k] == "-":
                            ss[k] = "T"
        for i in range(2, R - 2):
            if ss[i] == "-" and connected(i - 2, i + 2):
                u = CA[i] - CA[i - 2]
                v = CA[i + 2] - CA[i]
                cosang = (u @ v) / (np.linalg.norm(u) * np.linalg.norm(v))
                if math.degrees(math.acos(min(1, max(-1, cosang)))) > 70.0:
                    ss[i] = "S"
        out[f] = ss
        if progress is not None:
            progress((f + 1) / len(frames))
    return out


DSSP_NAMES = {"H": "alpha helix", "G": "3-10 helix", "I": "pi helix", "E": "strand", "B": "bridge",
              "T": "turn", "S": "bend", "-": "coil"}


# =============================================================================
#  8. SASA
# =============================================================================
def _sphere_points(n):
    """Golden-spiral points on the unit sphere."""
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    theta = np.pi * (1 + 5 ** 0.5) * k
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)


def sasa(coords, radii, probe=1.4, n_points=100, box=None):
    """Solvent accessible surface area per atom (A^2), Shrake-Rupley."""
    X = np.asarray(coords, dtype=float)
    R = np.asarray(radii, dtype=float) + probe
    unit = _sphere_points(n_points)
    rmax = R.max()
    i, j, _d = pairs_within(X, X, 2 * rmax, box)
    keep = i != j
    i, j = i[keep], j[keep]
    order = np.argsort(i, kind="stable")
    i, j = i[order], j[order]
    starts = np.searchsorted(i, np.arange(len(X)), side="left")
    ends = np.searchsorted(i, np.arange(len(X)), side="right")
    area = np.zeros(len(X))
    for a in range(len(X)):
        nbrs = j[starts[a]:ends[a]]
        pts = X[a] + R[a] * unit
        if len(nbrs):
            d = pts[:, None, :] - X[nbrs][None, :, :]
            if box is not None:
                d -= box * np.round(d / box)
            buried = ((d * d).sum(axis=2) < (R[nbrs] ** 2)[None, :]).any(axis=1)
            free = n_points - int(buried.sum())
        else:
            free = n_points
        area[a] = 4.0 * math.pi * R[a] ** 2 * free / n_points
    return area


# =============================================================================
#  9. Mean square displacement and self-diffusion
# =============================================================================
def unwrap(frames, boxes):
    """Removes the jumps of atoms crossing the periodic boundaries."""
    X = np.asarray(frames, dtype=float).copy()
    if boxes is None:
        return X
    for f in range(1, len(X)):
        box = np.asarray(boxes[f], dtype=float)
        d = X[f] - X[f - 1]
        X[f] = X[f - 1] + d - box * np.round(d / box)
    return X


def mean_square_displacement(frames, max_lag=None, boxes=None, remove_com=True):
    """MSD(lag) averaged over atoms and time origins (A^2)."""
    X = unwrap(frames, boxes)
    if remove_com:
        X = X - X.mean(axis=1, keepdims=True)
    F = len(X)
    max_lag = F - 1 if max_lag is None else min(int(max_lag), F - 1)
    lags = np.arange(max_lag + 1)
    msd = np.zeros(max_lag + 1)
    for lag in lags[1:]:
        d = X[lag:] - X[:-lag]
        msd[lag] = (d * d).sum(axis=2).mean()
    return lags, msd


def diffusion_coefficient(times_ps, msd, fit_from=0.2, fit_to=0.8):
    """D (1e-5 cm^2/s) from the slope of MSD between the given fractions of the
    time range (Einstein: MSD = 6 D t). Returns (D, slope, intercept, i0, i1)."""
    n = len(times_ps)
    i0, i1 = int(fit_from * (n - 1)), max(int(fit_to * (n - 1)), int(fit_from * (n - 1)) + 2)
    t, y = np.asarray(times_ps[i0:i1 + 1]), np.asarray(msd[i0:i1 + 1])
    slope, intercept = np.polyfit(t, y, 1)
    D = slope / 6.0 * 1e-4 / 1e-5        # . A^2/ps -> cm^2/s (1 A^2/ps = 1e-4 cm^2/s), in 1e-5 cm^2/s
    return D, slope, intercept, i0, i1


# =============================================================================
#  10. Interaction energy
# =============================================================================
def lj_pair_function(form, eps, sig):
    """Per-pair LJ parameters for a combination form: returns a function of
    (i_atoms, j_atoms, r) giving the energy (kJ/mol)."""
    eps = np.asarray(eps, dtype=float)
    sig = np.asarray(sig, dtype=float)
    form = (form or "OPLS").upper()

    def energy(i, j, r):
        e = np.sqrt(eps[i] * eps[j])
        if form == "AMBER":
            rmin = 0.5 * (sig[i] + sig[j])
            s6 = (rmin / r) ** 6
            return e * (s6 * s6 - 2.0 * s6)
        s = np.sqrt(sig[i] * sig[j])
        s6 = (s / r) ** 6
        return 4.0 * e * (s6 * s6 - s6)
    return energy


def interaction_energy(frames, group_a, group_b, charges, lj_energy=None, excluded=None, scaled14=None,
                       scale_elec14=0.5, lj14_energy=None, boxes=None, cutoff=None, b_groups=None, progress=None):
    """Coulomb and LJ energy (kJ/mol) between group A and group B per frame.

    excluded: set of (i, j) pairs (i < j; 1-2, 1-3) left out; scaled14: set of
    1-4 pairs, with Coulomb scaled by scale_elec14 and LJ from lj14_energy
    (the 1-4 parameters of the force field). cutoff=None: all pairs.
    b_groups: array mapping every atom of B to a group (e.g. residue) for the
    decomposition. Returns (elec (F,), lj (F,), per-group mean total {g: kJ/mol}).
    """
    A = np.asarray(group_a, dtype=int)
    B = np.asarray(group_b, dtype=int)
    q = np.asarray(charges, dtype=float)
    excluded = excluded or set()
    scaled14 = scaled14 or set()
    elec, ljs = [], []
    per_group = {}
    for f, X in enumerate(frames):
        box = None if boxes is None else boxes[f]
        if cutoff is None and box is None:
            ia = np.repeat(np.arange(len(A)), len(B))
            jb = np.tile(np.arange(len(B)), len(A))
            d = X[A][ia] - X[B][jb]
            r = np.sqrt((d * d).sum(axis=1))
        else:
            ia, jb, r = pairs_within(X[A], X[B], cutoff or 1e6, box)
        i, j = A[ia], B[jb]
        fe = np.ones(len(i))
        fl = np.ones(len(i))
        is14 = np.zeros(len(i), dtype=bool)
        if excluded or scaled14:
            for k in range(len(i)):
                key = (i[k], j[k]) if i[k] < j[k] else (j[k], i[k])
                if key in excluded:
                    fe[k] = fl[k] = 0.0
                elif key in scaled14:
                    fe[k], fl[k], is14[k] = scale_elec14, 0.0, True
        r = np.maximum(r, 0.5)
        e_el = COULOMB_KJ * q[i] * q[j] / r * fe
        e_lj = (lj_energy(i, j, r) * fl) if lj_energy is not None else np.zeros(len(i))
        if is14.any() and lj14_energy is not None:
            e_lj[is14] += lj14_energy(i[is14], j[is14], r[is14])
        elec.append(e_el.sum())
        ljs.append(e_lj.sum())
        if progress is not None:
            progress((f + 1) / len(frames))
        if b_groups is not None:
            g = np.asarray(b_groups)[jb]
            tot = np.bincount(g, weights=e_el + e_lj)
            for gi in np.nonzero(tot)[0]:
                per_group[int(gi)] = per_group.get(int(gi), 0.0) + tot[gi]
    n = max(len(frames), 1)
    return np.array(elec), np.array(ljs), {g: v / n for g, v in per_group.items()}


# =============================================================================
#  11. Dihedrals and rotamers
# =============================================================================
def dihedral(frames, i, j, k, l):
    X = np.asarray(frames, dtype=float)
    b0 = X[:, i] - X[:, j]
    b1 = X[:, k] - X[:, j]
    b2 = X[:, l] - X[:, k]
    b1n = b1 / np.linalg.norm(b1, axis=1, keepdims=True)
    v = b0 - (b0 * b1n).sum(axis=1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(axis=1, keepdims=True) * b1n
    x = (v * w).sum(axis=1)
    y = (np.cross(b1n, v) * w).sum(axis=1)
    return np.degrees(np.arctan2(y, x))


# . standard side-chain chi angles (IUPAC atom names)
CHI_ATOMS = {
    "ARG": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")],
    "ASN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "ASP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "CYS": [("N", "CA", "CB", "SG")],
    "GLN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")],
    "GLU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")],
    "HIS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")],
    "ILE": [("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")],
    "LEU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "LYS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")],
    "MET": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD")],
    "PHE": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "PRO": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")],
    "SER": [("N", "CA", "CB", "OG")],
    "THR": [("N", "CA", "CB", "OG1")],
    "TRP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "TYR": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "VAL": [("N", "CA", "CB", "CG1")],
}
for _alias, _name in (("HID", "HIS"), ("HIE", "HIS"), ("HIP", "HIS"), ("HSD", "HIS"), ("HSE", "HIS"),
                      ("HSP", "HIS"), ("CYX", "CYS"), ("ASH", "ASP"), ("GLH", "GLU"), ("LYN", "LYS")):
    CHI_ATOMS[_alias] = CHI_ATOMS[_name]


def rotamer_state(angles):
    """g+ (0, 120), t (120, 240) i.e. around 180, g- (240, 360), per IUPAC sign convention."""
    a = np.mod(np.asarray(angles), 360.0)
    return np.where(a < 120.0, "g+", np.where(a < 240.0, "t", "g-"))


def rotamer_summary(angles):
    states = rotamer_state(angles)
    pops = {s: float((states == s).mean()) for s in ("g+", "t", "g-")}
    transitions = int((states[1:] != states[:-1]).sum())
    return pops, transitions


# =============================================================================
#  12. Statistics
# =============================================================================
def block_average(series):
    """Flyvbjerg-Petersen blocking analysis of a correlated series.

    Returns (levels, sem): levels is a list of (n_blocks, sem, sem_error) for
    successive halvings; sem is the value at the first plateau (two successive
    levels agreeing within the error bar), or the largest SEM if no plateau is
    reached (the series is too short - a conservative estimate).
    """
    x = np.asarray(series, dtype=float)
    levels = []
    while len(x) >= 4:
        n = len(x)
        sem = math.sqrt(x.var(ddof=1) / n)
        levels.append((n, sem, sem / math.sqrt(2.0 * (n - 1))))
        x = 0.5 * (x[0:2 * (n // 2):2] + x[1:2 * (n // 2):2])
    if not levels:
        return levels, float("nan")
    for k in range(len(levels) - 1):
        if levels[k + 1][1] - levels[k][1] < levels[k + 1][2]:
            return levels, levels[k + 1][1]
    return levels, max(l[1] for l in levels)


def autocorrelation(series, max_lag=None):
    """Normalised autocorrelation function and the integrated correlation time
    (in frames, summed until the ACF first drops below zero)."""
    x = np.asarray(series, dtype=float)
    x = x - x.mean()
    n = len(x)
    max_lag = n // 2 if max_lag is None else min(max_lag, n - 1)
    var = (x * x).mean()
    if var == 0:
        return np.ones(1), 0.0
    nfft = 1 << (2 * n - 1).bit_length()
    fx = np.fft.rfft(x, nfft)
    acf = np.fft.irfft(fx * np.conj(fx), nfft)[:max_lag + 1] / (n - np.arange(max_lag + 1)) / var
    tau = 0.5
    for v in acf[1:]:
        if v <= 0:
            break
        tau += v
    return acf, tau


def describe(series):
    x = np.asarray(series, dtype=float)
    _levels, sem = block_average(x)
    _acf, tau = autocorrelation(x)
    return {"n": len(x), "mean": float(x.mean()), "std": float(x.std(ddof=1)) if len(x) > 1 else 0.0,
            "min": float(x.min()), "max": float(x.max()), "sem_block": float(sem), "tau_frames": float(tau),
            "n_independent": float(len(x) / max(2 * tau, 1.0))}


# =============================================================================
#  B-factors and average structure
# =============================================================================
def fluctuations(aligned_frames):
    """RMSF (A), isotropic B (A^2, 8 pi^2/3 <dr^2>) and anisotropic U (n, 3, 3)
    of superposed frames, as pDynamo's CoordinateFluctuations."""
    X = np.asarray(aligned_frames, dtype=float)
    mean = X.mean(axis=0)
    D = X - mean
    U = np.einsum("fni,fnj->nij", D, D) / X.shape[0]
    msf = np.trace(U, axis1=1, axis2=2)
    return {"average": mean, "rmsf": np.sqrt(msf), "b_factors": 8.0 * math.pi ** 2 / 3.0 * msf, "U": U}


# =============================================================================
#  Quasi-harmonic analysis
# =============================================================================
def quasi_harmonic(aligned_frames, masses, temperature=300.0):
    """Quasi-harmonic frequencies (cm^-1) and the vibrational entropy
    (J/(mol K), quantum harmonic oscillators, Andricioaei & Karplus 2001) from
    the mass-weighted covariance; also Schlitter's (1993) upper bound.

    The frames must be superposed: the six overall rotation/translation
    directions then have (near) zero variance and are dropped.
    """
    p = pca(aligned_frames, masses)
    lam = p["eigenvalues"] * AMU_KG * 1e-20        # . kg m^2
    lam = lam[lam > 1e-30 * AMU_KG]
    kT = KB_JK * temperature
    omega = np.sqrt(kT / lam)                       # . rad/s
    nu_cm = omega / (2.0 * math.pi * C_CM_S)
    x = HBAR_JS * omega / kT
    s_modes = x / np.expm1(x) - np.log1p(-np.exp(-x))
    S_qh = KB_JK * AVOGADRO * s_modes.sum()
    # . Schlitter: S < k/2 ln det(1 + kT e^2/hbar^2 M sigma)
    e2 = math.e ** 2
    S_sch = 0.5 * KB_JK * AVOGADRO * np.log1p(kT * e2 / HBAR_JS ** 2 * lam).sum()
    return {"frequencies_cm": nu_cm, "entropy_qh": float(S_qh), "entropy_schlitter": float(S_sch),
            "temperature": temperature, "modes": len(lam), "per_mode": KB_JK * AVOGADRO * s_modes}
