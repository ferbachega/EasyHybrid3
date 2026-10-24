---
slug: ts_search
title: Transition state search (Baker)
section: User guide
order: 6
kind: Manual
summary: The Baker saddle-point search, with a worked cyclohexane example.
---
# Transition state search (Baker) in EasyHybrid: user manual

> Draft, 2026-10-19. It describes the current development version of EasyHybrid, not 3.1.3.

## 1. What the tool does

**Simulate → Search Transition State (Baker)...** refines a structure that is already close to a
transition state into a **first-order saddle point**: a stationary point whose Hessian has exactly one
negative eigenvalue.

It runs pDynamo's `BakerSaddleOptimize_SystemGeometry`, the eigenvector-following method of
J. Baker (*J. Comput. Chem.* **7**, 385–395, 1986). It is the same function used in Example 11 of the
pDynamo book.

The method works in three steps:
1. A numerical Hessian is computed once, at the starting structure.
2. At every step the search goes **uphill** along one Hessian mode (the reaction coordinate) and
   **downhill** along all the others.
3. The Hessian is not recomputed. After each step it is updated from the change in the gradient.

The method **refines** a transition state; it does not find a reaction path. The starting structure
must already be near the saddle point. Good starting points are:
- the energy maximum of a PES scan (**Simulate → Reaction Coordinate Scans**);
- the highest image of an NEB path (**Simulate → Nudged Elastic Band**) or a CPR path (**Simulate → Conjugate Peak Refinement...**);
- a hand-built guess, with the forming and breaking bonds at intermediate lengths.

The search runs as an ordinary job in the Process Manager, with the same log, trajectory and import
handling as a geometry optimization.

---

## 2. The window

| Field | Default | Meaning |
|---|---|---|
| **Coordinates** | last set | Starting structure, chosen from the coordinate sets of the active system. |
| **Log frequency** | 1 | Write one line to the log every N steps: energy, RMS gradient and number of negative Hessian eigenvalues. |
| **Maximum iterations** | 250 | Maximum number of steps. If the search reaches it, the job is reported as **Finished (not converged)**. |
| **RMS Grad. Tolerance** | 0.001 | The search stops when the RMS gradient falls below this value (kJ mol⁻¹ Å⁻¹). |
| **Follow mode** | 0 | Hessian mode followed uphill, counted from the lowest eigenvalue (0 = lowest). |
| **Hessian updating** | BOFILL | How the Hessian is updated between steps: BOFILL (recommended) or POWELL. |
| **Recompute Hessian every** | 0 | Recompute the exact Hessian every N steps. 0 means it is computed once at the start and then only updated. |
| **Save trajectory** | off | Writes the steps to a `.ptGeo` folder (or `.dcd`, `.mdcrd`, `.xyz`) together with `output.log`. |

Every field has a tooltip with the same explanation.

### Follow mode

At the first step the search climbs the chosen Hessian mode. After that it keeps following the same
mode by its overlap with the previous step, even if the eigenvalue order changes.

- **0** (default) is right when the starting structure is close to the transition state. There the
  lowest mode, with negative curvature, is the reaction coordinate.
- **1, 2, ...** are for cases where the reaction coordinate is a higher mode, for example when starting
  from a minimum where the lowest modes are methyl rotations. Run **Simulate → Normal Modes** on the
  starting structure and pick the mode that moves the atoms of the reaction.

### Recompute Hessian every

By default the Hessian is computed once, at the start, and then only updated from the gradients. With
N > 0 the exact Hessian is recomputed every N steps. This needs fewer steps and is more robust far from
the transition state, but each recalculation is expensive.

pDynamo has no analytic Hessian, so the exact Hessian comes from finite differences, as in Normal
Modes. Each recalculation costs 6 gradient evaluations per free atom. That is cheap with MM or
semi-empirical QC and expensive with ab initio or DFT.

For cyclohexane (section 3), the search takes 123 steps with N = 0, 86 with N = 10, 80 with N = 5 and
27 with N = 1. All four reach the same transition state.

### Two pDynamo3 problems fixed in EasyHybrid

EasyHybrid runs its own copy of pDynamo's Baker optimizer (`src/pdynamo/p_methods/baker_saddle.py`),
because pDynamo3's version has two problems in these options:
- with Follow mode > 0 the search stops after 3 steps with an internal error
  (`AbsoluteMaximumIndex`);
- the "analytic Hessian" steps compute nothing, because pDynamo3 has no analytic Hessian for molecular
  systems, so they were wasted.

The algorithm is otherwise unchanged.

### Why only BOFILL and POWELL

We ran the cyclohexane example (section 3) with each update scheme:

| Hessian updating | Result |
|---|---|
| **BOFILL** | Transition state in 123 steps, one imaginary frequency (−311 cm⁻¹). |
| **POWELL** | Same transition state in 135 steps. |
| BFGS | Did not converge in 250 steps. BFGS keeps the Hessian positive, which is wrong for a saddle point. |
| Murtagh–Sargent | "Converged" back to the chair, a minimum, not a transition state. |

BFGS and Murtagh–Sargent are therefore not offered in the window. BFGS also fails with the Hessian
recomputed every 5 steps: after 250 steps it is back at the chair.

---

## 3. Worked example: cyclohexane chair → half-chair transition state

This is Example 11 of the pDynamo book, run from the EasyHybrid interface. The search starts from the
chair conformation, a minimum, and climbs to the saddle point of the ring flip. The calculation takes a
few seconds.

**Files used** (both ship with pDynamo3):
- structure: `$PDYNAMO3_HOME/examples/book/data/mol/cyclohexane_chair.mol`
- OPLS parameters: `$PDYNAMO3_HOME/parameters/forceFields/opls/bookSmallExamples`

### Step 1. Load the system

1. Click **New Project** on the toolbar.
2. Set the force field type to **OPLS**.
3. In the OPLS parameter folder chooser, select the `bookSmallExamples` folder. The default is `protein`,
   which has no cyclohexane types.
4. Choose `cyclohexane_chair.mol` as the coordinate file and import.

The energy of the chair is **82.56 kJ/mol** (**Simulate → Single Point Energy**).

### Step 2. Run the search

1. Open **Simulate → Search Transition State (Baker)...**.
2. In **Coordinates**, keep the imported chair structure.
3. Keep all defaults: 250 iterations, tolerance 0.001, Follow mode 0, BOFILL.
4. Optionally tick **Save trajectory** to keep the steps.
5. Click **Run**.

### Step 3. Check the result

When the job finishes, the "Finished!" dialog opens and the new structure is added as a coordinate set.

| Quantity | Expected value |
|---|---|
| Steps | 123 |
| Final RMS gradient | < 0.001 |
| Energy | ≈ 177.8 kJ/mol, i.e. **+95.3 kJ/mol** above the chair |
| Last line of the optimizer output in the log | `Minimization converged.` |

pDynamo writes "Minimization converged." for saddle-point searches too; it means the convergence
criterion was met.

Then confirm that the structure is a first-order saddle point:
1. Run **Simulate → Normal Modes** on the new coordinate set.
2. The log must list **exactly one negative (imaginary) frequency**, about **−311 cm⁻¹**.

The energy barrier here is a property of the small OPLS test parameters used in the pDynamo book. It is
not a reference value for the real cyclohexane ring flip.

---

## 4. When the search does not converge

If the search reaches **Maximum iterations** before the gradient falls below the tolerance:
- the Process Manager shows **Finished (not converged)**;
- the "Finished!" dialog shows a warning with the number of steps and the final RMS gradient;
- the same warning is written to the log.

**The structure from such a run is not a transition state.** For example, cyclohexane stopped after 50
steps with tolerance 0.1 has an RMS gradient of 3.1 and two imaginary frequencies.

What to do:
- **Restart from the last structure** with more iterations. The search keeps no state between runs, so
  a restart recomputes the Hessian at the new point. This often helps.
- **Start from a better guess**, for example the maximum of a finer PES scan.
- **Recompute the Hessian** every few steps (**Recompute Hessian every** = 1 to 10).
- **Try POWELL** instead of BOFILL.

If pDynamo itself reports an error, the job is marked **Error!** and no structure is imported.

---

## 5. After the transition state: the reaction path

From a converged transition state, **Simulate → IRC / Reaction Path...** follows the imaginary mode
downhill in both directions, towards reactants and products. It checks again that the Hessian has
exactly one negative eigenvalue, so a structure from a non-converged search will be rejected there.

The Reaction Path tool has no convergence test of its own. Look at the whole path, not only at its
final structure.

---

## 6. Practical notes

- **Cost.** Each numerical Hessian (the starting one, and every recalculation) takes 6N gradient
  calculations for N free atoms. That is cheap
  with MM or semi-empirical QC, but expensive with DFT on large systems. For QC/MM systems, keep the
  number of free atoms small (fix distant atoms) so the Hessian stays small.
- **Fixed atoms** are honoured. The search runs only over the free coordinates.
- **Restraints** are part of the energy. Remove PES-scan restraints before searching, or the result is
  a saddle point of the restrained surface.
- **Units.** Energies are in kJ/mol and the gradient tolerance is in kJ mol⁻¹ Å⁻¹, as everywhere in
  pDynamo.

## 7. Reference

J. Baker, "An algorithm for the location of transition states", *J. Comput. Chem.* **7**, 385–395
(1986).
