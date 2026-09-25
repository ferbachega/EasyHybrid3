# Translating a DYFF (UFF) Molecular-Mechanics Model into CHARMM-Format Topology and Parameter Files for NAMD

## 1. Overview and rationale

pDynamo3's `MMModelDYFF` implements DYFF, an in-house reformulation of Rappé
and co-workers' Universal Force Field (UFF) [1], extended with metal-organic
framework parameters from Addicoat et al. [2] and Coupry et al. [3]. DYFF
reproduces UFF's original bonded functional forms — including its
non-harmonic, Fourier-in-cosine valence-angle and out-of-plane (inversion)
potentials — essentially unmodified.

NAMD [4] has no native reader for pDynamo3's own MM representation. Its only
generic force-field input path accepts either an AMBER `prmtop`/`inpcrd` pair
or a CHARMM-format topology/parameter/coordinate triple (`.psf`/`.prm`/`.pdb`),
both of which restrict every bonded term to a fixed, small set of
closed-form functional expressions (harmonic bonds and angles, a truncated
Fourier series for proper dihedrals, and a harmonic improper term). Because
those functional forms are not, in general, capable of representing DYFF's
own potential exactly, running a DYFF-parameterized system in NAMD requires a
translation step, not merely a file-format conversion.

This document describes that translation in full: which terms carry over
exactly, which require a controlled approximation, how the approximation is
constructed and validated, and the specific numerical/software pitfalls
encountered while producing files NAMD accepts without silently discarding
or misassigning parameters. The translation is implemented as a
self-contained Python module, `src/util/dyff_to_charmm.py`, operating
directly on pDynamo3's in-memory `System`/`MMModelDYFF` objects (as restored
from a YAML- or pickle-serialized pDynamo3 System) — not on any intermediate
text export — so that every numeric value used below is read from the same
data structures pDynamo3 itself uses to evaluate the system's energy and
forces, rather than re-derived from tabulated UFF constants.

The guiding design principle throughout is: **carry over exactly whatever
can be carried over exactly, and make every remaining approximation
explicit, quantified, and reproducible from first principles**, rather than
silently refitting the whole force field or approximating terms that in
fact need no approximation at all.

## 2. Functional forms

### 2.1 DYFF

Every bonded term container in pDynamo3's MM state exposes its parameters as
one of two closed forms, confirmed directly against the compiled energy
evaluators (`CosineTermEnergies.c`, `HarmonicBondContainer.c`) rather than
assumed from the literature:

- **Bonds** (`HarmonicBondContainer`): a plain harmonic,
  $$E_{\text{bond}}(b) = k_b\,(b-b_0)^2$$
  with no leading factor of one-half — confirmed from the C evaluator's
  `energy += fc * disp * disp` accumulation.

- **Angles, proper dihedrals, and out-of-plane (inversion) terms**
  (`CosineAngleContainer`, `CosineDihedralContainer`,
  `CosineOutOfPlaneContainer`): a general truncated Fourier series in the
  cosine of the relevant angle,
  $$E(x) = \sum_{p} c_p \cos(p\,x)$$
  where $x$ is the valence angle $\theta$, the proper dihedral angle
  $\phi$, or DYFF's own out-of-plane angle $\omega$ (defined in §2.3 below),
  and each term type's $\{(c_p,p)\}$ set is looked up per interacting
  atom-type tuple. This is confirmed directly from the persisted container
  state (`__getstate__` returns exactly these `(coefficient, period)` pairs)
  and from the compiled evaluator, which computes $E$, $dE/dx$, and the
  gradient by explicit summation over these terms (`CosineTermEnergy_Angle`,
  `_Dihedral`, `_OutOfPlane`).

For angles specifically, this generalizes UFF's own case-by-case
trigonometric expressions (single-term for linear centres, a two-term
$\propto[1-\cos(3\theta)]$ form for trigonal-planar centres, and a
three-term form fitted to the equilibrium angle and a reference force
constant for the general case) into one uniform representation, which is
exploited below (§3.2) to locate each type's true equilibrium angle and
local curvature purely numerically, without special-casing any of UFF's own
hybridization rules.

### 2.2 CHARMM (as read by NAMD)

$$
\begin{aligned}
E_{\text{bond}} &= K_b (b-b_0)^2 \\
E_{\text{angle}} &= K_\theta (\theta-\theta_0)^2 \\
E_{\text{dihedral}} &= K_\chi \left[\,1+\cos(n\chi-\delta)\,\right] \\
E_{\text{improper}} &= K_\psi (\psi-\psi_0)^2 \\
E_{\text{vdW}}(r_{ij}) &= \varepsilon_{ij}\!\left[\left(\frac{R_{\min,ij}}{r_{ij}}\right)^{12} - 2\left(\frac{R_{\min,ij}}{r_{ij}}\right)^{6}\right]
\end{aligned}
$$

with CHARMM's conventional combining rule $\varepsilon_{ij}=\sqrt{\varepsilon_i\varepsilon_j}$
(geometric) and $R_{\min,ij}=R_{\min,i}+R_{\min,j}$, i.e. the arithmetic
mean of the two per-type $R_{\min}/2$ values, doubled.

By construction, none of `angle`, `improper`, or `vdW` (for a heteroatomic
pair) can represent DYFF's corresponding term exactly; `bond` and
`dihedral` can, as shown next.

### 2.3 DYFF's out-of-plane angle vs. CHARMM's improper dihedral

For an out-of-plane centre $j$ bonded to $i$, $k$, $l$, DYFF's angle
$\omega$ is the angle between the $j\!\to\! i$ bond vector and the unit
normal $\hat{n}$ of the plane spanned by the $j\!\to\! k$ and $j\!\to\! l$
bond vectors:

$$\hat n = \widehat{(\mathbf r_k-\mathbf r_j)\times(\mathbf r_l-\mathbf r_j)}\,,\qquad
\cos\omega = \hat n\cdot\widehat{(\mathbf r_i-\mathbf r_j)}$$

CHARMM's improper "dihedral" $\psi$, by contrast, is the ordinary four-body
torsion angle of the same four atoms taken as a chain $i\!-\!j\!-\!k\!-\!l$
— algebraically identical to a proper-dihedral calculation, only its
physical role differs. These are genuinely different angular coordinates of
the same four atoms (§4.4).

## 3. Exact and approximate term-by-term mapping

| Term | Mapping | Nature |
|---|---|---|
| Bond stretch | direct | **exact** |
| Proper dihedral | algebraic reduction | **exact**, conditional on a property verified against the real data (§3.3) |
| Valence angle bend | curvature-matched harmonic | approximation (equilibrium value exact; shape away from it is not) |
| Out-of-plane / improper | curvature-matched harmonic + numeric angle-coordinate conversion | approximation (per-instance geometric correction) |
| van der Waals | direct self-pair; explicit `NBFIX` for every heteroatomic pair | **exact**, via `NBFIX` (§3.5) |
| Partial charges / 1–4 electrostatic scaling | direct | **exact** |

### 3.1 Bond stretching

Because DYFF's and CHARMM's bond potentials share the identical functional
form $k(b-b_0)^2$, the mapping is a direct unit conversion with no
re-fitting:

$$K_b\ [\text{kcal mol}^{-1}\text{Å}^{-2}] = \frac{k_b\ [\text{kJ mol}^{-1}\text{Å}^{-2}]}{4.184},\qquad b_0\ [\text{Å}]\ \text{unchanged.}$$

### 3.2 Valence angle bending

For each DYFF angle-parameter type, the true equilibrium angle $\theta_0$
is located as the global minimum of $E(\theta)=\sum_p c_p\cos(p\theta)$ over
$\theta\in(0,\pi)$, by a fine grid scan (guarding against the multiple
local extrema a general $N$-term series can have) followed by Newton
refinement on $dE/d\theta=0$. Because $\theta_0$ is a true minimum,
$dE/d\theta|_{\theta_0}=0$, and the Taylor expansion about it is

$$E(\theta_0+\delta) = E(\theta_0) + \tfrac12 E''(\theta_0)\,\delta^2 + O(\delta^4)\,,\qquad
E''(\theta)= -\sum_p p^2 c_p\cos(p\theta)\,.$$

Matching the quadratic coefficient to CHARMM's $K_\theta(\theta-\theta_0)^2$
fixes

$$K_\theta = \frac{E''(\theta_0)}{2}\quad(\text{kJ mol}^{-1}\text{rad}^{-2}\to\text{kcal mol}^{-1}\text{rad}^{-2}\text{ by }\div 4.184)\,,$$

which is the **unique** harmonic term reproducing DYFF's local stiffness
(hence the same small-amplitude bending frequency) at the exact DYFF
equilibrium angle. It necessarily departs from the true, generally
anharmonic and multi-well DYFF potential away from $\theta_0$; §5 quantifies
the resulting total-angle-energy discrepancy for the validation system
(+4.7%, evaluated at an unequilibrated, far-from-minimum starting
structure — see discussion there on why this is closer to an upper bound
than a typical operating value).

No Urey–Bradley term is added ($K_{UB}=0$): DYFF has no corresponding 1–3
term, so introducing one would add curvature CHARMM's own angle-only
representation cannot attribute correctly.

### 3.3 Proper dihedrals

UFF's dihedral term for a bond $j$–$k$ has the closed form (Rappé et al.
[1], eq. 16)

$$E(\phi) = \tfrac12 V_{jk}\left[\,1-\cos(n\phi_0)\cos(n\phi)\,\right]
= \underbrace{\tfrac12 V_{jk}}_{c_0} \;+\; \underbrace{\left(-\tfrac12 V_{jk}\cos(n\phi_0)\right)}_{c_1}\cos(n\phi)\,,$$

i.e. exactly a two-term DYFF Fourier series with $c_1=\pm c_0$, the sign
fixed by $\cos(n\phi_0)$. UFF's own parameterization always chooses $n$ and
$\phi_0$ such that $\cos(n\phi_0)=\pm1$ exactly (an integer multiple of
$\pi$), which was **verified programmatically against every dihedral
parameter type actually present in the target system** (not assumed): every
one is a 2-term series with $|c_0|=|c_1|$ to machine precision. The
translation code enforces this check and raises rather than silently
mis-converting a type for which it should ever fail (e.g. for a future,
different DYFF-parameterized system).

Given $|c_0|=|c_1|$, the elementary identity $\cos(x-\pi)=-\cos x$ (which
holds unconditionally, independent of $n$) gives an **exact** algebraic
match to CHARMM's dihedral form:

$$
E(\phi)=c_0\bigl(1-\cos n\phi\bigr)=c_0\bigl(1+\cos(n\phi-\pi)\bigr)\ \Rightarrow\ K_\chi=|c_0|,\ \delta=180^\circ \quad(c_1=-c_0)
$$
$$
E(\phi)=c_0\bigl(1+\cos n\phi\bigr) \Rightarrow K_\chi=|c_0|,\ \delta=0^\circ \quad(c_1=+c_0)
$$

with $K_\chi\,[\text{kcal mol}^{-1}] = |c_0|\,[\text{kJ mol}^{-1}]/4.184$.
This is not an approximation: the two representations are the same
function.

### 3.4 Out-of-plane terms / CHARMM impropers

DYFF's out-of-plane energy, like its angle term, is a general series
$E(\omega)=\sum_p c_p\cos(p\omega)$; its true minimum $\omega_0$ and local
curvature $E''(\omega_0)$ are obtained by the identical numeric procedure
described in §3.2 (restricted to $\omega\in(0,\pi/2)$, the physically
distinct range given the term's built-in $[0,\pi/2]$ periodicity by
construction — see the DYFF source's own out-of-plane reformulation note).

Because $\omega$ and CHARMM's improper $\psi$ are different angular
coordinates of the same four atoms (§2.3), curvature-matching cannot be
done directly in $\omega$; it must be re-expressed in $\psi$. Both
coordinates reach their extremal value at exactly the same configuration —
the planar one, $\omega=\pi/2$, $\psi\in\{0,\pi\}$ — because four coplanar
points force the two spanning planes in each angle's own definition to
coincide, an elementary geometric fact independent of either force field.
Writing $\omega=\omega(\psi)$ locally near that shared extremum and using
$dE/d\omega|_{\omega_0}=0$, the chain rule gives

$$\left.\frac{d^2E}{d\psi^2}\right|_{\text{eq}} = \left.\frac{d^2E}{d\omega^2}\right|_{\omega_0}\left(\left.\frac{d\omega}{d\psi}\right|_{\text{eq}}\right)^{2}\,,$$

so

$$K_\psi = \frac{1}{2}\left.\frac{d^2E}{d\omega^2}\right|_{\omega_0}\left(\left.\frac{d\omega}{d\psi}\right|_{\text{eq}}\right)^{2}\,.$$

The purely geometric factor $d\omega/d\psi$ depends on the local bond
geometry around the specific improper centre and is therefore evaluated
**numerically, per real improper instance**, not assumed constant: atom $i$
is displaced by $\pm h$ along the current DYFF normal $\hat n$, and both
$\omega$ and $\psi$ are recomputed at $\mathbf r_i\pm h\hat n$; their
central-difference ratio gives $d\omega/d\psi$ at that instance's actual
(possibly slightly non-planar, since the target structure is not
pre-equilibrated) geometry. Because chemically equivalent instances of one
DYFF out-of-plane parameter type can still sit at slightly different local
geometries in a real structure, the ratio is averaged over every instance
sharing that type before combining it with the type's curvature, yielding
one $K_\psi$ per type — matching CHARMM's static, per-type parameter model.
$\psi_0$ is assigned per instance as whichever of $0^\circ/180^\circ$ the
instance's actual geometry is closer to (both are equally valid planar
solutions; which applies is a matter of the arbitrary atom-ordering
convention used to list the four atoms, not of the underlying physics).

### 3.5 Nonbonded (van der Waals) interactions

`MMModelDYFF`'s Lennard-Jones parameters are stored in `LJForm.OPLS` form:
combining is **geometric** for *both* $\varepsilon$ and $\sigma$,
$\varepsilon_{ij}=\sqrt{\varepsilon_i\varepsilon_j}$,
$\sigma_{ij}=\sqrt{\sigma_i\sigma_j}$ — matching UFF's own combining rule
[1] exactly. NAMD/CHARMM's built-in combining rule uses the same geometric
rule for $\varepsilon$ but an **arithmetic** one for $R_{\min}$
($R_{\min,ij}=R_{\min,i}+R_{\min,j}$), which would silently give the wrong
pairwise potential for any heteroatomic pair if left to its own default.

This is resolved exactly, not approximately: every unordered pair of atom
types actually present in the system is written as an explicit `NBFIX`
override, computed directly with DYFF's own (geometric) combining rule
rather than relying on NAMD's default:

$$\varepsilon_{ij}=\sqrt{\varepsilon_i\varepsilon_j}\,,\qquad
R_{\min,ij} = 2^{1/6}\sqrt{\sigma_i\sigma_j}\,,$$

each computed twice — once for the regular nonbonded interaction and once
for the separately-scaled 1–4 interaction, since DYFF's own 1–4
Lennard-Jones table (`ljParameters14`) uses a distinct, halved
$\varepsilon$ at the same $\sigma$ (see §3.6). Self-pairs ($i=i$) require no
`NBFIX`, since the per-type `NONBONDED` line already gives the exact value.
$R_{\min}$ throughout is derived from $\sigma$ via the standard 12–6
relation $R_{\min}=2^{1/6}\sigma$ (the separation at the potential minimum),
not assumed equal to $\sigma$ itself.

### 3.6 Partial charges and 1–4 scaling

Atomic partial charges (CM5, in this case) are copied through unmodified.
`MMModelDYFF` exposes its own 1–4 scaling factors explicitly —
`electrostaticScale14` and `lennardJonesScale14` — which for the validation
system are both $0.5$, **read directly from the model rather than assumed**
from a generic UFF convention (the two need not coincide, and did not need
to be guessed). The electrostatic value maps directly onto NAMD's
`1-4scaling` keyword; the Lennard-Jones value is already realized concretely
via the separate 1–4 `NONBONDED`/`NBFIX` columns (§3.5), which is how
CHARMM/NAMD apply 1–4 vdW scaling in general (as fixed, explicit
per-type/per-pair values rather than a single global multiplier).

Exclusions require no special handling beyond NAMD's standard
`exclude scaled1-4` policy: pDynamo3's own exclusion pair list was verified
(by direct count) to equal the union of 1–2, 1–3, and 1–4 pairs implied by
the system's bond connectivity, and its 1–4 interaction list to equal
exactly the set of dihedral end-atom pairs — i.e. the standard organic
connectivity-derived exclusion pattern `scaled1-4` already assumes, with no
DYFF-specific exclusions to reproduce separately.

## 4. Atom-type refinement

DYFF's own parameter selection is not a pure function of the two (or three,
or four) interacting atoms' coarse hybridization labels alone: for bonded
terms, it also depends on the discrete **bond order** of the specific
bond(s) spanned by the term, via UFF's bond-order-dependent natural bond
length expression (a logarithmic correction term, Rappé et al. [1], eq. 3).
CHARMM/NAMD's parameter lookup has no equivalent run-time concept — it is
purely a function of each atom's single, static type as recorded once in
the `.psf` file.

This mismatch was **detected empirically**, not anticipated in the abstract:
mapping each atom directly onto its coarse DYFF hybridization label (e.g.
`C:Tri`, `O:Tri`) and grouping the resulting bond/angle instances by
interacting-type tuple revealed that two or more numerically distinct
parameter sets collapsed onto the same CHARMM type key for several tuples
in the validation system — for example, an aromatic-ring C–C bond
($b_0=1.385$ Å) and an exocyclic ring-to-carbonyl C–C bond ($b_0=1.464$ Å)
both map to the coarse pair (`C:Tri`, `C:Tri`); an ester C–O single bond
($b_0=1.346$ Å) and a carbonyl C=O double bond ($b_0=1.220$ Å) both map to
(`C:Tri`, `O:Tri`). Left unresolved, CHARMM's parameter lookup can only ever
see one of the two colliding values for a given type pair, silently
discarding the other for whichever fraction of bonds/angles happens to use
it — in this system, roughly 14.5% of bonds and 23% of angles by instance
count.

The resolution is mechanical rather than chemistry-specific: every atom's
coarse type label is refined with a binary suffix, `_D`/`_S`, recording
whether *any* of its real covalent bonds has a DYFF-computed bond order
$\geq 1.9$ (i.e. is a genuine double bond), using DYFF's own bond-order
function (`DYFFBondOrder`) rather than an independently reimplemented rule.
This refinement was verified, by exhaustive re-checking of every bond,
angle, dihedral, and improper parameter-type group, to eliminate every
numeric collision in the validation system; it is a no-op (produces an
unused type variant) for atom types that cannot structurally carry a double
bond (e.g. DYFF's own tetrahedral carbon/oxygen types, or hydrogen),
so it generalizes safely to systems where no such collision exists at all.

One further, related subtlety follows directly from this refinement and
was likewise caught empirically: a single underlying DYFF parameter *type*
index can, after refinement, correspond to *more than one* distinct literal
atom-type quadruple/triple/pair — an angle or dihedral type shared by
several chemically equivalent placements will fan out into different
`_D`/`_S` label combinations depending on which specific atom occupies
which position in a given instance. NAMD's parameter file must contain a
line for *every distinct literal combination actually used in the
topology*, not merely one representative combination per underlying DYFF
type; the implementation groups parameter output by the observed literal
label tuple, not by the DYFF type index, for this reason (confirmed the
hard way: an earlier version of the translation, grouping by DYFF type
index alone, produced a `.prm` file NAMD rejected outright with "UNABLE TO
FIND DIHEDRAL PARAMETERS").

## 5. Implementation and validation

### 5.1 File format notes

Topology (`.psf`) and parameter (`.prm`) files are written directly by the
translation code (not via any third-party topology tool). The `.psf` uses
CHARMM's "extended"/X-PLOR convention of string (rather than numeric-index)
atom types, matching the type labels used in the accompanying `.prm`.
Because NAMD's PSF reader tokenizes each `ATOM` record on whitespace, every
field is written with an explicit separating space regardless of its
individual width; an earlier version that instead relied on fixed-width,
unseparated columns silently merged adjacent fields whenever a value (e.g.
a five-character sanitized DYFF type label) filled its nominal column
exactly, corrupting the following field. Coordinates are written as a
standard, strict-column PDB (name/altLoc/resName/chain/resSeq in their
canonical columns), since NAMD associates `.psf` and `.pdb` records purely
by atom count and order, not by name matching, but strict columns keep the
file usable by other common tools (e.g. VMD).

### 5.2 Validation against the source model

For a representative system (7200 atoms; a hydrated lipid-bilayer-like
assembly with two ester/aromatic-containing lipid residues and TIP3-like
water, 6200 bonds / 8600 angles / 10800 dihedrals / 2100 out-of-plane
terms, 5 coarse / 7 refined atom types), the translated files were loaded
into NAMD 3.0.1 and their reported bonded energies at the starting
(unrelaxed) geometry were compared against pDynamo3's own energy
decomposition for the identical coordinates:

| Term | pDynamo3 (kcal mol⁻¹) | NAMD (kcal mol⁻¹) | Difference |
|---|---:|---:|---:|
| Bond | 1830.3702 | 1830.3826 | +0.00 % |
| Dihedral | 270.6504 | 270.6508 | +0.00 % |
| Improper | 0.0071 | 0.0071 | +0.36 % |
| Angle | 4866.1596 | 5093.9927 | +4.68 % |

The bond, dihedral, and improper agreement (all $<0.4\%$, essentially at
floating-point/rounding precision) confirms the exact-mapping and
per-instance geometric-conversion procedures (§3.1, 3.3, 3.4) are correctly
implemented end-to-end, not merely correct in derivation. The angle
discrepancy is consistent with, and of the expected order for, the
curvature-matched harmonic approximation (§3.2) evaluated far from
equilibrium — the test geometry is deliberately the system's raw,
unrelaxed starting structure (chosen to also stress-test the nonbonded
terms, see below), so this figure should be read as closer to an upper
bound on the approximation's practical impact than a typical value once a
structure is minimized/equilibrated, at which point most angles sit much
closer to their own $\theta_0$ and the quadratic approximation is
correspondingly tighter.

The same starting geometry carries severe steric overlaps unrelated to the
translation itself (independently confirmed via pDynamo3's own nonbonded
energy, $\sim\!5\times10^{10}$ kJ mol⁻¹), which is expected of a freshly
assembled, unminimized structure. This was exploited as a further
functional test: a 5000-step energy minimization (PME electrostatics,
orthorhombic $40\times40\times50$ Å periodic cell) using the translated
files converged smoothly, with **no warnings or errors**, from a total
potential energy of $+1.05\times10^{7}$ kcal mol⁻¹ (dominated by the
initial steric overlap) to $-1.53\times10^{4}$ kcal mol⁻¹, in 4.9 s of
wall-clock time on 4 CPU threads. Continuing from that minimized structure,
a 10000-step (10 ps) isothermal–isobaric (NPT, Langevin thermostat/piston,
298 K, 1.01325 bar) production run likewise completed with no
warnings or errors (94.4 ns day⁻¹ on the same hardware), with the
instantaneous temperature settling to $293$–$298$ K and the periodic cell
volume fluctuating within $\sim\!1$–$4\%$ of its initial value, i.e. stable
NPT behaviour rather than divergence.

## 6. Scope and limitations

- The angle and out-of-plane/improper approximations (§3.2, §3.4) reproduce
  DYFF's equilibrium geometry and local (harmonic) stiffness exactly, but
  not its full anharmonic, periodic shape; they should be expected to
  under- or over-estimate energies for large-amplitude bending/inversion
  excursions relative to the native DYFF potential. No Urey–Bradley
  coupling is introduced to compensate, for the reason given in §3.2.
- The atom-type refinement (§4) is a bond-order-only disambiguation. It
  was verified sufficient for the bonded-parameter set of the validation
  system; a system whose DYFF parameterization depends on some other
  discrete structural feature not reflected in bond order alone would
  require extending the same principle (refine by whatever discrete,
  mechanically-derivable feature DYFF's own parameter assignment actually
  used) rather than assuming coverage.
- `NBFIX` coverage is exhaustive over the atom types *actually present* in
  a given system, not a general-purpose DYFF/UFF parameter table; a
  different DYFF-parameterized system regenerates its own complete,
  self-consistent set from scratch (the translation is run per-system, not
  from a pre-built library).
- PSF `NGRP` (charge-group) records are written as a single, degenerate
  group spanning the whole system. NAMD's own nonbonded evaluation is
  atom-based, not charge-group-based, so this has no effect on any reported
  energy or force; it is noted here only because it departs from what a
  CHARMM-lineage topology tool (e.g. `psfgen`) would normally emit.

## References

[1] Rappé, A. K.; Casewit, C. J.; Colwell, K. S.; Goddard, W. A., III;
    Skiff, W. M. "UFF, a Full Periodic Table Force Field for Molecular
    Mechanics and Molecular Dynamics Simulations." *J. Am. Chem. Soc.*
    **1992**, *114*, 10024–10035.

[2] Addicoat, M. A.; Vankova, N.; Akter, I. F.; Heine, T. "Extension of the
    Universal Force Field to Metal-Organic Frameworks." *J. Chem. Theory
    Comput.* **2014**, *10*, 880–891.

[3] Coupry, D. E.; Addicoat, M. A.; Heine, T. "Extension of the Universal
    Force Field for Metal-Organic Frameworks." *J. Chem. Theory Comput.*
    **2016**, *12*, 5215–5225.

[4] Phillips, J. C. et al. "Scalable molecular dynamics on CPU and GPU
    architectures with NAMD." *J. Chem. Phys.* **2020**, *153*, 044130.

## Implementation

The full translation is implemented in `src/util/dyff_to_charmm.py`
(EasyHybrid3 repository), a self-contained Python module with no
EasyHybrid-specific dependencies beyond pDynamo3 itself, usable either as a
library (`convert(yaml_path, output_prefix)`) or as a command-line tool.
