"""
charge_blocks
=============

Integer charge blocks of an MM system (the Charge Block Inspector, System >
Charges). No GTK.

    blocks = residue_blocks(system)                  # residues and their atoms
    rows   = analyze(charges, blocks, tolerance)     # charge, target, deviation
    new    = correct(charges, rows, method)          # corrected charges

The charge of every residue is brought to an integer target (the nearest
integer, or one set by the user) and the difference is spread over its atoms:

  even          the same shift on every atom: dq_i = d / n
  proportional  proportional to |q_i|: dq_i = d |q_i| / sum |q|. Atoms without
                charge (often nonpolar H, dummy atoms) keep zero; falls back
                to even when the residue has no charge at all.
"""

import numpy as np

METHODS = [("even", "Same shift on every atom"),
           ("proportional", "Proportional to |q| (zero-charge atoms kept)")]


def residue_blocks(system):
    """ [{'key', 'chain', 'label', 'name', 'indices'}] from the pDynamo sequence (one block
        with every atom when the system has no sequence). """
    blocks = []
    sequence = getattr(system, "sequence", None)
    if sequence is not None and getattr(sequence, "children", None):
        for entity in sequence.children:
            for component in entity.children:
                indices = [atom.index for atom in component.children]
                if not indices:
                    continue
                label = component.label
                blocks.append({"key": len(blocks), "chain": entity.label, "label": label,
                               "name": label.split(".")[0], "indices": indices})
    if not blocks:
        blocks.append({"key": 0, "chain": "", "label": "all atoms", "name": "all",
                       "indices": list(range(len(system.atoms)))})
    return blocks


def analyze(charges, blocks, tolerance=1.0e-4, targets=None):
    """ One row per block: charge, target (nearest integer or targets[key]), deviation
        (charge - target) and whether it needs a correction (|deviation| > tolerance). """
    q = np.asarray(charges, dtype=float)
    targets = targets or {}
    rows = []
    for b in blocks:
        total = float(q[b["indices"]].sum())
        target = float(targets.get(b["key"], round(total)))
        deviation = total - target
        rows.append(dict(b, charge=total, target=target, deviation=deviation, n=len(b["indices"]),
                         needs=abs(deviation) > tolerance, user_target=b["key"] in targets))
    return rows


def correct(charges, rows, method="even"):
    """ New charges with every block that needs it at its target. """
    q = np.asarray(charges, dtype=float).copy()
    for r in rows:
        if not r["needs"]:
            continue
        idx = np.asarray(r["indices"], dtype=int)
        d = r["deviation"]
        if method == "proportional":
            weights = np.abs(q[idx])
            total = weights.sum()
            if total > 1.0e-12:
                q[idx] -= d * weights / total
                continue
        q[idx] -= d / len(idx)
    return q


def summary(charges, new_charges, rows):
    q0 = np.asarray(charges, dtype=float)
    q1 = np.asarray(new_charges, dtype=float)
    dq = q1 - q0
    changed = [r for r in rows if r["needs"]]
    worst = max(changed, key=lambda r: abs(r["deviation"])) if changed else None
    big_target = [r for r in changed if abs(r["target"] - round(r["charge"])) > 0.5]
    return {"total_before": float(q0.sum()), "total_after": float(q1.sum()),
            "n_blocks": len(rows), "n_changed": len(changed),
            "n_atoms_changed": int(np.sum(np.abs(dq) > 1.0e-12)),
            "max_atom_shift": float(np.abs(dq).max()) if len(dq) else 0.0,
            "worst": worst, "user_targets": big_target}
