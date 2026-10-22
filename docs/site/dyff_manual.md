# DYFF in EasyHybrid: user manual

> Hidden page (draft, 2026-10-16; updated 2026-10-21). It describes the current development
> version of EasyHybrid, not 3.1.3.

## 1. What DYFF is

**DYFF**, the *Dynamo Force Field*, is pDynamo's general force field. It comes from the
Universal Force Field (UFF) of Rappé et al. and includes the metal-organic framework extensions
of Addicoat et al. and Coupry et al.

DYFF has two useful properties:
- **Atom types come from the chemistry alone:** element, hybridization and bond orders. Any
  molecule, including ones with metals, can be typed without a residue template.
- **It has no electrostatics of its own.** The charges must come from elsewhere.

Plain UFF is meant for rough geometries, not for production simulations. EasyHybrid adds
several layers on top of pDynamo's `dyff-1.0` parameters to make DYFF usable in production
simulations. We call the result **modified DYFF**.

| Layer | What it adds |
|---|---|
| **Charges** | CM5 charges from GFN1-xTB for small molecules; charges from the residue library for proteins. |
| **Polar-hydrogen LJ** | H bonded to N, O or S gets a small Lennard-Jones type (σ = 1.0 Å). Hydrogen bonds then do not sit inside the repulsive wall of plain UFF hydrogens. |
| **xTB refinement** | Bond and angle force constants from an xTB Hessian (modified Seminario method), calibrated against DFT for each chemical environment. |
| **DFT torsion library** | Torsion profiles for each chemical environment, fitted to ~18 000 B3LYP-D3BJ/DZVP torsion scans. Rotatable bonds not in the library get a relaxed xTB scan in the Builder (Path A), and optionally in *Prepare DYFF System* (Path B, section 4.1). |
| **Residue library** | Ready-made, force-matched parameters and charges for the standard amino acids and their protonation states. |
| **Metal sites** | Metal–donor bonds, DYFF metal types, and xTB parameters for the metal and its ligands. |
| **Water and ions** | TIP3P charges and oxygen LJ for water; formal charges for ions. |

The result is an ordinary pDynamo MM model. Minimization, molecular dynamics, scans, umbrella
sampling and QC/MM all work with it.

---

## 2. Which path to follow

| You have | Follow |
|---|---|
| A small molecule or an isolated ligand | **Path A**: Builder and DYFF Parameters |
| A protein, with or without ligands, cofactors, solvent or metal ions | **Path B**: *System › Force Field › Define Force Field › Prepare DYFF System…* |
| An isolated metal complex | Path A |
| A metalloprotein | Path B |

**Requirement:** xTB must be configured in EasyHybrid (*Preferences*). RDKit is no longer
required.

---

## 3. Path A: small molecule

**A1. Open the molecule in the Builder.** Toolbar button *Builder*, then:
- **New** to draw, or to build from the fragment and structure libraries;
- **Edit existing** to work on a copy of a system already loaded (mol2, sdf, pdb…). Bond orders
  are perceived automatically. The original system is only replaced when you accept at the
  end.

**A2. Check the chemistry.** This is the most important step: DYFF types, environment codes and
library matches all depend on it.
- **Bond orders.** Check aromatic rings, carboxylates and charged groups.
- **Formal charges.** Use *Protonate* / *Deprotonate*.
- **Hydrogens.** Use *Add Hydrogens…*, or the *Adjust hydrogen count* switch.

**[CONFIRM]** The *DYFF Parameters* button was removed from the Builder toolbar on 2026-10-21.
Steps A3–A6 need a new entry point (or a rewrite) before this page is published.

**A3. Open *DYFF Parameters…*** from the Builder sidebar. The window lists:
- the DYFF type of each atom;
- every bond, angle, dihedral and out-of-plane term, colored by its strain in the current
  geometry (green, yellow, red).

Clicking a row highlights the atoms.

**A4. Charges.** Click *CM5 charges…* (GFN1-xTB) to give the model electrostatics. Polar
hydrogens get their small LJ type automatically. *Clear charges* removes them.

- **Use charges** for solvent, intermolecular interactions and complexes.
- **Skip them** for isolated conformations in the gas phase: charges make gas-phase torsion
  profiles slightly worse.
- **Do A4 before A5.** The torsion table used by the refinement depends on whether the model
  has charges.

**A5. *Refine with xTB (Hessian + torsions)…*** The dialog offers:
- **Hessian** (recommended): GFN2-xTB Hessian at the xTB minimum; bond and angle constants
  calibrated to DFT.
- **Torsions from the DFT library** (default on): library profiles where the chemical
  environment is known, relaxed xTB scans (24 points, run in parallel) for the rest.
- **Charge** and **multiplicity**.
- **ALPB implicit water**: switched on automatically when the molecule has a formal charge.
- **Whole molecule** or **Selected atoms only**:
  - the selection is refined as a capped sub-system with a neighborhood of 3 bonds by default;
  - only terms whose atoms are all selected change.

The whole refinement is one *Undo* step. The status bar reports:
- how many terms were calibrated;
- how many torsions came from the library and how many from scans;
- any warning.

**A6. Optional fine tuning.**
- Use the strain colors and the *attention* filter to find problem terms.
- On a single dihedral, *Fit with xTB scan…* plots the xTB and DYFF profiles and fits that
  dihedral.
- *Clear xTB refinement* goes back to the DYFF rules.

**A7. Close the Builder** and answer **Yes** to *Incorporate with DYFF and customizations?*.
The pDynamo system is created, or replaced in *Edit existing*. It carries:
- the DYFF model and the term edits;
- the charges and the polar-H LJ;
- a non-bonded model.

Answering **No** gives plain DYFF.

---

## 4. Path B: protein or metalloprotein

**B1. Open the structure** with *File › Import Data…*, or open the PDB directly.
- A PDB without CONECT records is read without bonds; the preparation creates the chemistry.
- For NMR ensembles, only the first model is used.

**B2. *System › Force Field › Define Force Field › Prepare DYFF System…*** and choose
**DYFF modified** (the other two choices are plain DYFF rules without xTB). A system without
bonds opens the preparation window below; a system that already has bonds and bond orders
opens the shorter *DYFF + residue library* dialog, with the same parameter options. The
structure is analyzed for:
- residues and inferred protonation;
- disulfides and chain breaks;
- metal sites (metal plus donors, found by distance);
- ligands, truncated residues, and covalently attached groups.

**B3. Choices in the window**

- **Covalently attached residues** (peptide caps ACE/NME/NH2, covalent ligands): removed by
  default; the chain is left with a free terminus.
- **Truncated residues:** missing atoms are placed with standard geometry, choosing the
  rotamer that avoids clashes. Check them and minimize.
- **Metal sites:** edit the model, charge and multiplicity for each site.
  - **Bonded** model: metal–donor bonds are created.
    - Coordinating Cys becomes thiolate.
    - His is protonated on the non-coordinating nitrogen.
    - Asp and Glu are deprotonated.
    - Coordinated water is left unbonded.
  - **Ionic** model: the metal is a point charge with LJ. This is the default for Mg and Ca.
  - A bonded site with fewer than 2 donors in the file becomes ionic, with a warning. This
    happens when the site is completed by crystal symmetry.
- **Ligands:** protonation at the chosen **pH**.
- **Method** for out-of-library groups: GFN2-xTB or g-xTB. Use GFN2-xTB with ALPB for charged
  sites, because g-xTB has no implicit solvent.
- **CM5 scale** and the **system name**.
- **Scan with xTB and fit the flexible bonds not in the torsion library** (off by default):
  see section 4.1.

Typical multiplicities: Zn(II) 1, Cu(II) 2, high-spin Fe(II) 5, high-spin Fe(III) 6.

**B4–B5. Automatic build.** A new system is created with:
- hydrogens, bond orders and residue names;
- metal–donor bonds and DYFF metal types (for example `Zn:Tet:+2`, `Ca:Oct:+2`);
- **standard residues** from the force-matched residue library;
- **out-of-library groups**, each refined in a capped sub-system: the metal with its donor
  residues, ligands, and unusual residues. Each group gets:
  - CM5 charges, using the site multiplicity;
  - an xTB Hessian, calibrated to DFT where an environment code exists;
  - torsions from the DFT library;
- **water and ions:** TIP3P charges and oxygen LJ for water, formal charges for ions.

Expect 1–5 min per metal site with GFN2-xTB; g-xTB is 2–3× slower.

**Alternative:** *System › Force Field › Define Force Field › Prepare OPLS System…*, then
*Prepare DYFF System…* › DYFF modified. This route does not accept Zn, Cu, Fe and similar ions.

### 4.1 Optional: automatic xTB torsion scans of flexible bonds

By default, the out-of-library groups take their torsions **only from the DFT torsion
library**. A rotatable bond whose chemical environment is not in the library keeps the generic
DYFF torsion rule, and nothing else is done with it. The final message reports the coverage:
"*X of Y rotatable bonds from the DFT torsion library*"; the remaining Y − X bonds kept the rule.

Tick **Scan with xTB and fit the flexible bonds not in the torsion library** to give those
bonds their own torsion:

1. **Which bonds.** Only bonds of the atoms outside the residue library (ligands, termini,
   modified residues, metal sites), and only **flexible** ones. Skipped as rigid or not worth a
   scan:
   - bonds in rings;
   - double, triple and aromatic bonds;
   - amide-like bonds, C(=O)–N, C(=S)–N, C(=N)–N with a planar nitrogen (amides, ureas,
     carbamates, amidines);
   - symmetric terminal rotors: CH₃, NH₂, NH₃⁺;
   - bonds already found in the DFT torsion library.
2. **The scan.** Each bond gets a relaxed GFN2-xTB scan (24 points, every 15°, with the same
   implicit-solvent setting as the charges), on the same capped sub-system used for the
   Hessian. The scans run in parallel.
3. **The fit.** The torsion around that bond is fitted, self-consistently, to the xTB profile,
   on top of the Seminario bonds and angles and the CM5 charges just computed (the same method
   as the Builder's *Refine with xTB*).
4. **Acceptance.** A fit is kept only if it lowers the error against the xTB profile by at
   least 10 %; otherwise the bond keeps the DYFF rule and the message says so. Fits that are kept
   but still have an error above 2 kcal/mol are reported with a warning.

The final message adds: "*xTB torsion scans: N of M flexible bond(s) fitted (RMSE vs xTB, mean:
a → b kcal/mol)*", plus the warnings. The choice is stored with the system, so re-applying the
library after an edit in the Builder repeats the scans.

**Cost.** One scan costs from seconds to a few minutes, depending on the size of the sub-system.
In our tests: ATP, 4 bonds in 160 s in total (Hessian included); folic acid, 1 bond in 67 s.

**What to expect.** The fits always improved on the DYFF rule in our tests, but not always by
much:

| Molecule | Bonds scanned | RMSE vs xTB, rule → fit | xTB barrier |
|---|---|---|---|
| ATP | the 4 P–O bonds of the phosphate bridges | 6–14 → 4.5–8 kcal/mol | ~6–7 kcal/mol |
| Folic acid | CH₂–pterin | 3.7 → 2.8 kcal/mol | 3.5 kcal/mol |
| Adenosine, flavone | none: the library covers every rotatable bond | – | – |

A large remaining error means that the relaxed xTB profile is not a pure torsion: along the
scan, internal hydrogen bonds form and break (phosphates, polyols), or the rest of the molecule
relaxes in a way the MM model reproduces differently. A torsion alone cannot absorb that. Check
these bonds before production runs.

**When to use it.** Ligands and modified residues with flexible bonds that matter for your
question (a linker, a side chain that changes conformation in the binding site), when the
coverage message shows bonds without library torsions. Leave it off for proteins with only
standard residues (nothing to scan) and when speed matters more.

---

## 5. Transition metals: what to expect

- DYFF types the metals, and the xTB Hessian gives the bonds and angles around them. In
  sub-systems, a metal always keeps its ligands and never gets a capping hydrogen.
- **Not validated against DFT:** the DFT libraries contain no metals, so these terms keep
  their xTB values. Compare with crystal structures before a production run.
- **Metal-containing cofactors (heme, Fe–S clusters) are not supported yet**; the window warns
  about them.
- **Weak axial bonds** are too short with GFN2-xTB (type-1 copper: Cu–S(Met) 2.8 Å in the
  crystal, 2.4 Å after preparation). Prefer g-xTB for copper.

---

## 6. After preparation

1. **Solvate** with *System › Solvate System…*: cube, rectangular box or sphere,
   neutralization, salt. The solvated system is a new system; the solute keeps its parameters.
   - Water gets TIP3P charges on DYFF bonded terms. This water model is **not validated**.
   - Equilibrate under NPT and check the density before a production run.
2. **Non-bonded model:** set automatically (full interactions without periodicity, cutoff in
   boxes). Change it with *System › Add › NB Model*.
3. **Simulate** with any EasyHybrid method.
4. **Save the project** (`.easy`): term edits, charges and the model are stored.

---

## 7. Decision guide

| Question | Recommendation |
|---|---|
| Charges or not? | Charges for solvent, intermolecular interactions, proteins and complexes. None for isolated gas-phase conformations. |
| Whole molecule or selection? | Whole molecule up to ~60 atoms. Above that, or when only one region matters (a ligand site, a modified residue), use the selection: with a 3-bond neighborhood it gave the same result as the whole molecule. |
| DFT library or xTB scans? | Keep the library on. It carries DFT data and costs nothing. In the Builder (Path A), scans fill the gaps automatically. In Path B the gaps keep the DYFF rule unless you tick the automatic scans (section 4.1). |
| Automatic xTB torsion scans (Path B)? | On for ligands or modified residues whose flexible bonds matter and are not in the library (see the coverage in the final message). Off for standard proteins and when speed matters. Check fits reported above 2 kcal/mol. |
| ALPB implicit water? | Whenever there is a formal charge. Isolated ions are poorly described by xTB in the gas phase. |
| DYFF or OPLS-AA for proteins? | DYFF + library: better forces and energies for distorted geometries, more protonation states, metals. OPLS-AA: better relative conformer energies and a validated water model. For long explicit-water production runs, OPLS-AA (or AMBER/CHARMM) is still the conservative choice. |
| Out-of-plane (improper) terms? | Leave them alone. Fitting them to the xTB Hessian made results worse, so they keep the DYFF rules. |

---

## 8. Troubleshooting

| Message or symptom | Cause and fix |
|---|---|
| "no bonds or no bond orders" | The system was read from a PDB without CONECT records. Prepare it first (Path B). |
| Few residues recognized | Non-standard residue names, missing heavy atoms, or unusual protonation. Unrecognized residues are treated as out-of-library groups. That works, but it is slower. |
| "N imaginary frequencies" | The geometry was not at an xTB minimum. Use *Use the xTB optimized geometry*, or optimize first. |
| "scan changed the bonding" | Proton transfer or bond breaking during a torsion scan. Use ALPB and check the protonation. That bond keeps the DYFF torsion. |
| "the fit does not improve on the DYFF rule …, rule kept" | Automatic scans (section 4.1): the fitted torsion was not at least 10 % better than the rule, so the rule was kept. |
| "fitted torsion(s) still above 2 kcal/mol RMSE" | Automatic scans (section 4.1): the xTB profile is not a pure torsion (internal hydrogen bonds?). The fit was kept because it is better than the rule; check those bonds. |
| "dihedral cut by the subsystem" | The four atoms of the scanned dihedral are not all inside the capped sub-system. That bond keeps the DYFF torsion. |
| "term(s) kept the DYFF value" | A non-physical constant came out of the Hessian, which is common for gas-phase anions. This is normal; ALPB usually removes it. |
| "Per-atom … edits were NOT applied" | The molecule changed after the refinement (an atom or bond was added). Run A5 again. |
| Missing table with charges | `dyff_bond_library_cm5.json` is missing; the table without charges is used instead. |

---

## 9. References

- **UFF:** Rappé, A. K.; Casewit, C. J.; Colwell, K. S.; Goddard, W. A.; Skiff, W. M.
  *J. Am. Chem. Soc.* **1992**, *114*, 10024–10035.
- **UFF metal extensions:**
  - Addicoat, M. A.; Vankova, N.; Akter, I. F.; Heine, T. *J. Chem. Theory Comput.*
    **2014**, *10*, 880–891.
  - Coupry, D. E.; Addicoat, M. A.; Heine, T. *J. Chem. Theory Comput.* **2016**, *12*,
    5215–5225.
- **xTB methods:**
  - GFN1-xTB: Grimme, S.; Bannwarth, C.; Shushkov, P. *J. Chem. Theory Comput.* **2017**,
    *13*, 1989–2009.
  - GFN2-xTB: Bannwarth, C.; Ehlert, S.; Grimme, S. *J. Chem. Theory Comput.* **2019**,
    *15*, 1652–1671.
  - ALPB: Ehlert, S.; Stahn, M.; Spicher, S.; Grimme, S. *J. Chem. Theory Comput.* **2021**,
    *17*, 4250–4261.
- **CM5 charges:** Marenich, A. V.; Jackson, S. V.; Cramer, C. J.; Truhlar, D. G.
  *J. Chem. Theory Comput.* **2012**, *8*, 527–541.
- **Seminario method:**
  - Seminario, J. M. *Int. J. Quantum Chem.* **1996**, *60*, 1271–1277.
  - Modified version: Allen, A. E. A.; Payne, M. C.; Cole, D. J. *J. Chem. Theory Comput.*
    **2018**, *14*, 274–281.
- **TIP3P water:** Jorgensen, W. L.; Chandrasekhar, J.; Madura, J. D.; Impey, R. W.;
  Klein, M. L. *J. Chem. Phys.* **1983**, *79*, 926–935.

**[CONFIRM]** Before publishing:
- cite the source of the DFT torsion scans (the QCArchive/OpenFF datasets used);
- add the DYFF paper, once it exists.
