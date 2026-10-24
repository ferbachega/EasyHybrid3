"""
et_pathways
===========

Electron-transfer coupling pathways (Beratan, Betts & Onuchic, Science 252,
1285, 1991; Beratan, Onuchic, Winkler & Gray, Science 258, 1740, 1992).

The donor-acceptor coupling is the product of the decays of the steps of
the best path through the protein:

    covalent bond      eps_C = 0.6
    hydrogen bond      eps_H = eps_C^2 exp[-1.7 (R - 2.8)]    R = heavy-atom D...A distance
    through space      eps_S = eps_C/2 exp[-1.7 (R - 1.4)]    R = atom-atom distance

The best path maximizes the product, i.e. it is the shortest path with
weights -ln(eps) (Dijkstra, from every donor atom to any acceptor atom).
sigma_l = ln(product) / ln(eps_C) * 1.4 A is the equivalent covalent
tunnelling length. pSimulation.ElectronTransferPathways is unfinished
("do not use"), so the model is implemented here.
"""

import heapq
import math

import numpy as np

from util import traj_analysis as TA

EPS_C = 0.6
DECAY = 1.7


def build_graph(elements, xyz, bonds, space_cutoff=5.0, hbond_cutoff=3.5, box=None):
    """Adjacency {i: [(j, eps, kind)]} with covalent, hydrogen-bond and through-space links."""
    n = len(elements)
    el = [TA.element_of(e) for e in elements]
    adj = {i: [] for i in range(n)}
    bonded = set()
    for i, j in bonds:
        adj[i].append((j, EPS_C, "C"))
        adj[j].append((i, EPS_C, "C"))
        bonded.add((min(i, j), max(i, j)))
    # . hydrogen bonds: donor heavy atom ... acceptor (heavy-atom distance)
    donors, acceptors = TA.donors_and_acceptors(el, bonds)
    hb = set()
    if donors and acceptors:
        D = np.array(sorted({d for d, _h in donors}))
        A = np.array(acceptors)
        ii, jj, rr = TA.pairs_within(xyz[D], xyz[A], hbond_cutoff, box)
        for a, b, r in zip(D[ii], A[jj], rr):
            if a == b or (min(a, b), max(a, b)) in bonded:
                continue
            eps = EPS_C ** 2 * math.exp(-DECAY * (r - 2.8))
            eps = min(eps, EPS_C)
            adj[int(a)].append((int(b), eps, "H"))
            adj[int(b)].append((int(a), eps, "H"))
            hb.add((min(a, b), max(a, b)))
    # . through space, between atoms not covalently bonded
    ii, jj, rr = TA.pairs_within(xyz, xyz, space_cutoff, box)
    for a, b, r in zip(ii, jj, rr):
        if a >= b:
            continue
        key = (int(a), int(b))
        if key in bonded or key in hb:
            continue
        eps = 0.5 * EPS_C * math.exp(-DECAY * (r - 1.4))
        adj[int(a)].append((int(b), eps, "S"))
        adj[int(b)].append((int(a), eps, "S"))
    return adj


def best_path(adj, donors, acceptors):
    """(coupling, [atoms], [kinds]) of the best path from any donor to any acceptor."""
    targets = set(acceptors)
    dist = {}
    prev = {}
    heap = [(0.0, d) for d in donors]
    for d in donors:
        dist[d] = 0.0
        prev[d] = (None, None)
    heapq.heapify(heap)
    found = None
    while heap:
        w, i = heapq.heappop(heap)
        if w > dist.get(i, math.inf):
            continue
        if i in targets:
            found = i
            break
        for j, eps, kind in adj[i]:
            nw = w - math.log(eps)
            if nw < dist.get(j, math.inf):
                dist[j] = nw
                prev[j] = (i, kind)
                heapq.heappush(heap, (nw, j))
    if found is None:
        return 0.0, [], []
    atoms, kinds = [found], []
    i = found
    while prev[i][0] is not None:
        j, kind = prev[i]
        atoms.append(j)
        kinds.append(kind)
        i = j
    atoms.reverse()
    kinds.reverse()
    return math.exp(-dist[found]), atoms, kinds


def sigma_length(coupling):
    """Equivalent covalent tunnelling length (A)."""
    if coupling <= 0:
        return float("inf")
    return math.log(coupling) / math.log(EPS_C) * 1.4


def pathways(frames, elements, bonds, donors, acceptors, boxes=None, space_cutoff=5.0, progress=None):
    """Best path in every frame. Returns a list of (coupling, atoms, kinds)."""
    out = []
    for f, X in enumerate(frames):
        box = None if boxes is None else boxes[f]
        adj = build_graph(elements, np.asarray(X, dtype=float), bonds, space_cutoff, box=box)
        out.append(best_path(adj, donors, acceptors))
        if progress is not None:
            progress((f + 1) / len(frames))
    return out
