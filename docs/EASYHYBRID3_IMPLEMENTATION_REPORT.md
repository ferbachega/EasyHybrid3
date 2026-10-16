# EasyHybrid3: A Graphical, Extensible Front End for QM/MM and Molecular Modelling with pDynamo3

**Internal technical report — implementation summary and detailed changelog**

---

## Abstract

EasyHybrid3 is a Python/GTK3 desktop application that provides an interactive graphical
front end to [pDynamo3](https://github.com/pdynamo/pDynamo3), a library for hybrid
quantum mechanics/molecular mechanics (QM/MM) simulation. The application couples
pDynamo3's simulation engine to a custom OpenGL molecular viewer, *vismol* (bundled as a
nested submodule), and adds a large body of original tooling on top of both: an
interactive molecule builder, docking front ends for AutoDock Vina and AutoDock-GPU,
enhanced-sampling and reaction-coordinate workflows (potential-energy surface scans,
umbrella sampling, WHAM, nudged elastic band, conjugate peak refinement), trajectory
and structural analysis (cartoon representations, radial distribution functions, RMSF,
Ramachandran plots, molecular surfaces, infrared spectra), and session/project
management facilities. This report documents the implemented functionality, the design
decisions behind it, and the validation methodology used throughout development. It is
organised in two parts: a narrative technical report (Sections 1–6) describing the
software as a whole, and a detailed, module-by-module changelog (Section 7) recording
individual features and bug fixes as they were implemented.

---

## Part I — Technical Report

### 1. Introduction

Hybrid QM/MM methods allow a small, chemically active region of a macromolecular
system to be described quantum-mechanically while the surrounding environment is
treated with a classical force field, making it possible to study reaction mechanisms
in realistic condensed-phase or enzymatic contexts at a fraction of the cost of a full
QM treatment. pDynamo3 provides the underlying algorithms (force fields, QM engines,
QM/MM coupling, geometry optimisation, reaction-path and free-energy methods) as a
Python library, but exposes no graphical interface of its own. EasyHybrid3 fills that
gap: it wraps pDynamo3's session/system model, adds a real-time 3D molecular viewer,
and layers a substantial set of interactive tools on top, so that system preparation,
QM/MM partitioning, simulation setup, and results analysis can all be carried out
without leaving a single graphical environment.

The project is under active, iterative development driven directly by hands-on usage:
new capabilities are typically motivated by a concrete modelling task the user is
carrying out, implemented, exercised against real molecular systems, and then refined
based on further use. This report reflects that history — most sections below describe
not only what a feature does today, but the concrete problem that motivated it and, in
several cases, real defects that were found and fixed by exercising the feature against
real molecules rather than synthetic test cases alone.

### 2. Software Architecture

EasyHybrid3 is organised around three cooperating layers:

- **pDynamo3 session layer** (`src/pdynamo/pDynamo2EasyHybrid/`) — wraps pDynamo3's
  `System` objects in a `pDynamoSession`, adding EasyHybrid-specific bookkeeping
  (`e_`-prefixed attributes: working folder, colour palette, restraint tables, QC/MM
  region tables, job history, and so on) alongside the numerical model itself.
- **Visualisation layer** (`src/graphics_engine/`, the *vismol* submodule) — an
  independent OpenGL/GTK molecular viewer, extended throughout this project with custom
  representations (sticks with true multi-bond rendering, ball-and-stick, a from-scratch
  cartoon ribbon renderer, molecular-surface meshes), a picking/selection system, and the
  low-level mouse/keyboard event pipeline the interactive Builder is built on.
  `VismolSession`/`VismolObject`/`Atom` form the visual data model; `VismolGLCore.render()`
  is the single per-frame draw/dispatch loop most interactive tools hook into.
- **Application layer** (`src/gui/`) — the GTK3 main window, menus, and the roughly
  seventy auxiliary windows (Glade-defined, one `.py` controller each) that implement
  individual tools: system setup, docking, sampling, analysis, and the Builder itself.
  `EasyHybridSession` (in `gui/eSession.py`) subclasses vismol's `VismolSession` and is
  the actual runtime session class, overriding several of its parent's methods.

A recurring architectural pattern is that a pDynamo `System` and its visual
`VismolObject` counterpart are linked by a shared integer key (`e_id`), with a system
capable of owning more than one visual object (for example, several docking poses, or a
derived molecular-surface mesh, sharing one `e_id`). A significant fraction of the work
described below — in particular the Builder's system-editing logic — is about keeping
these two parallel representations, and the auxiliary bookkeeping dictionaries that key
off them (treeview rows, per-system colour palettes, selection sets, undo stacks),
correctly synchronised as the user edits, clones, merges, prunes, or deletes systems and
objects.

### 3. Implementation

#### 3.1 Interactive Molecular Builder

The single largest body of work in the project is an interactive molecule builder,
allowing a user to construct or modify chemical structures directly in the 3D view
rather than only importing pre-built files. Two independent, code-duplicated
implementations were found early in the project's history; the dead one
(`src/graphics_engine/src/builder/`) was left untouched and all further work targeted
the live implementation (`src/gui/windows/builder/`).

**Core editing primitives.** Atom placement, deletion, bonding, and element/bond-order
changes are implemented as low-level operations in `atom_ops.py`
(`add_atom`/`remove_atom`/`add_bond`/`set_bond_order`/`set_atom_element`), each
correctly maintaining a dense, zero-based `atom_id` indexing scheme against which
positions, bonds, and every other per-atom array must stay aligned as atoms are added or
removed — a recurring source of subtle bugs (renumbering-after-removal hazards, all
resolved by re-reading each surviving atom's own `.atom_id` attribute, which
`remove_atom()` updates in place, rather than trusting a previously captured integer).
A dedicated undo stack captures full structural snapshots before each logical user
action (never at the level of individual low-level calls), so a single click-and-drag
gesture that internally issues several primitive operations still undoes as one step.

**Hydrogenation.** New or modified atoms have their hydrogen complement kept consistent
with standard valence via `adjust_hydrogens()`, which computes each atom's target
valence from its element and current formal charge, and both adds/removes hydrogens and
positions them using a VSEPR-style geometric rule shared with the heavy-atom skeleton
relaxation used by the "Clean Up" tool (a two-phase, non-force-field geometric
regularisation: bond lengths and angles first for heavy atoms, then hydrogens placed
against the now-relaxed skeleton). A per-session "Adjust hydrogen count" switch, wired
into every automatic post-edit hydrogenation call site (atom placement, bond dragging,
element replacement, bond-order changes, fragment attachment), allows this automatic
correction to be disabled, so that non-standard, deliberately over- or under-coordinated
species (for example a tetravalent, positively charged ammonium nitrogen) can be built
by hand, one bond at a time, without the tool silently reverting each addition.

**Bond creation and editing.** Bonds may be created by dragging from an existing atom
(with the new bond's order inferred from drag distance when no explicit order is
selected, and a live colour-coded preview ring), by a keyboard shortcut acting on two
picked atoms, or by a dedicated "Bond Order" click tool that changes an existing bond's
order in place. Bond order/aromaticity selection is exposed as a small set of real,
indicator-drawing radio buttons in the sidebar (an earlier "flat toggle button" style was
replaced at the user's request); an aromaticity option that existed briefly was later
removed in favour of the "Adjust hydrogen count" switch's own manual-construction
capability.

**Protonation state.** A `formal_charge` attribute was added to the visual `Atom` model,
kept deliberately distinct from the pre-existing `charge` attribute (which holds an MM
partial charge populated only after a force field is assigned, an entirely different
quantity). The "Atom Types" inspector window (below) was extended with an editable
integer "Charge" column and "Protonate"/"Deprotonate" actions that change an atom's
formal charge by one unit and let the existing hydrogen-adjustment machinery add or
remove the corresponding proton, positioned by the same VSEPR logic used everywhere
else. The formal charge is propagated into the real pDynamo `Atom.formalCharge` field
when the underlying pDynamo `System` is (re)built, so downstream force-field typing and
total-charge accounting see the correct, chemically meaningful state rather than a
purely cosmetic label.

**Fragment and structure libraries.** Two complementary mechanisms allow larger
chemical motifs to be introduced in one action. The *fragment library* holds small
functional groups (carboxylic acid, amide, isocyanate, methoxy, alkyl chains, phenyl,
furyl, and others) authored as ordinary, chemically complete `.mol2` files with exactly
one atom deliberately left one bond short of its standard valence; that single
under-valent atom is auto-detected as the attachment point, so authoring a new fragment
is a matter of drawing a complete molecule and deleting one terminal hydrogen. Attaching
a fragment replaces a chosen hydrogen atom with the fragment, oriented via a single
rigid rotation that aligns the fragment's own "missing substituent" direction (computed
by the same VSEPR machinery used for hydrogenation) with the real bond direction being
replaced. The *structure library* holds complete, standalone molecules (populated with
several hundred real compounds retrieved via the PubChem PUG REST API and processed with
Open Babel, organised into medicinal-chemistry scaffold classes) that can be placed as a
whole at an arbitrary point, with `.mol2` "ar" (aromatic) bond codes explicitly
re-Kekulised through the project's own bond-order perception routine rather than
collapsed to a uniform bond order. Both library windows offer a live, independently
rendered 3D preview.

Frequently used fragments are additionally exposed directly in the Builder sidebar as a
quick-pick button row; picking either an element or a quick-pick fragment shares a single
GTK radio-button group, so exactly one "what will the next click place" choice is active
at a time, and the sidebar's general "Add" tool is armed automatically by either kind of
pick.

**Editing an already-loaded system.** Beyond building from a blank canvas, an
already-loaded system (with one or more associated visual objects, for example several
superposed docking poses) can be edited directly: the user selects which visual object
to edit (if more than one is present), is warned if an already-assigned force field or
QM model would be discarded by the edit, and the Builder then operates on a disposable
clone rather than the original. When editing finishes, an atom-count-preserving edit is
folded back onto the original system in place (preserving its `e_id` and every other
piece of bookkeeping keyed against it, such as selections, restraints, and the QM/MM
region table); an edit that changes the atom count is instead kept as its own new,
independent system, since folding a changed atom count back in place would silently
desynchronise any of that index-keyed state. This "edit in place" / "edit into a new
system" split is a genuine correctness distinction, not a stylistic one, and — after an
earlier, more heavy-handed implementation was found to discard *other* visual objects
belonging to the same system regardless of which outcome the session eventually took —
was reworked to guarantee the original system, and every object in it, is left
completely unmodified whenever the "new system" outcome occurs. For visual clarity while
an edit is in progress, every other currently visible object (including the very object
being edited, since the Builder itself shows a working clone of it) is temporarily
hidden rather than removed, and made visible again, unconditionally, the instant the
Builder session ends.

**Auxiliary inspection and manipulation tools.** An "Atom Types" window shows, for the
current viewing selection, each atom's DYFF-perceived force-field type (computed live,
without requiring a force field to have actually been assigned yet) and its formal
charge, both editable by hand for cases the automatic perceiver gets wrong. A "Transform
Selection" window exposes six sliders (translation and rotation about X/Y/Z) acting on
the current selection about its own centroid. A dedicated "Dihedral Angle" window drives
live rotation of one side of a bond about a user-picked four-atom dihedral. A periodic
table dialog (covering the full main-group and transition-metal range, with UFF-derived
radii and electronegativity reused from the same source DYFF itself cites) supplements
the sidebar's quick-pick element buttons for less common elements. A one-click "Optimize
(DYFF)" action runs a real, bounded (100-step) DYFF force-field minimisation on a
disposable scratch copy of the current structure and copies the relaxed coordinates
back, without requiring the user to leave the Builder or manually assign a force field
to the working object.

**Viewing-selection interaction.** The application's general-purpose shift-click/
shift-drag "viewing selection" mechanism (used by several non-Builder tools to pick
atoms for analysis) is suppressed — both its visual rendering and the underlying
selection state itself — while the Builder's own "Editing: ON" mode is active, and
restored automatically the moment editing is paused or the Builder is closed, so the two
independent selection concepts (which atoms will the next Builder click affect, versus
which atoms are marked for some other tool's use) never visually or functionally
conflict with one another.

#### 3.2 Molecular Docking

Front ends for two docking engines were integrated directly into the session/treeview
model rather than as external, disconnected tools. An AutoDock Vina interface (ported
from an existing standalone tool) and an AutoDock-GPU interface (built by reusing the
Vina tool's own code, with grid-file preparation adjusted to avoid a real AutoGrid4
path-length limitation) both write PDBQT input, launch the docking engine as a tracked
background job, and import the resulting poses. A subsequent revision removed the
intermediate PDBQT round-trip for the common case, importing docking poses directly into
`VisMol` objects backed by a real, in-memory pDynamo `System` rather than temporary
files, with a mechanism to pin which trajectory frame a given pose-derived object
displays.

#### 3.3 Enhanced Sampling and Reaction-Coordinate Methods

A family of related tools compute reaction coordinates and drive the system along them.
A potential-energy-surface (PES) scan tool supports one- and two-dimensional scans over
distance, multi-distance, and dihedral reaction coordinates, with a "simple" mode
(fixed-shape reaction coordinates) and an "advanced" mode (an arbitrary weighted sum of
distances) later merged into a single window behind a mode toggle, sharing one
unmodified computational backend. The same merge pattern was applied to Umbrella
Sampling, which also gained genuine dihedral-restraint support (`RestraintDihedral`,
with correct 360° periodicity) rather than only distance-based restraints, and to an
"Energy Refinement" tool that measures (without scanning) a reaction coordinate value
along an already-computed trajectory. A WHAM (weighted histogram analysis method)
post-processing tool turns a set of umbrella-sampling windows into a potential-of-
mean-force profile, with an interactive statusbar readout and curve-identification
feature layered on top of the underlying plotting widgets, and automatic PNG export of
the resulting histograms and free-energy profile. Reaction-path methods are represented
by a nudged-elastic-band-style path tool and a conjugate peak refinement (CPR) tool built
by generalising the NEB tool's own initial-path construction to a new refinement
algorithm; both carry the caveat that the underlying pDynamo3 optimiser used to refine
the path has no built-in convergence check, so a returned path's final frame alone
cannot be assumed to represent a true saddle point without inspecting the accompanying
energy profile.

#### 3.4 Trajectory and Structural Analysis

A from-scratch cartoon (ribbon) representation was implemented independently of
pDynamo3/vismol's own prior attempts, including secondary-structure-dependent geometry,
correct handling of live trajectory playback, and a substantial (order-of-magnitude)
vectorised performance pass, with an optional caching mode trading memory for a further
roughly two-fold speed-up. A general molecular-surface tool computes van der Waals,
solvent-accessible, and solvent-excluded (Connolly) surfaces via a signed-distance-field
plus marching-cubes construction, coloured by the same per-atom colours the rest of the
viewer uses (including a mode that colours only carbon atoms by a system's custom
reference colour, leaving other elements on the standard CPK palette), with triangle-
area computation exposed for eventual partial-surface queries. Standard structural
analyses — radial distribution functions, per-residue RMSF, and Ramachandran plots — and
an infrared-spectrum tool (built on a corrected `QCModelXTB.DipoleMoment` implementation,
since the system-level dipole moment does not sum contributions from mixed MM/QC regions)
round out the analysis toolset.

#### 3.5 QM/MM Setup and External Engine Integration

Beyond pDynamo3's own native QM models, EasyHybrid3 provides setup front ends for
external engines: NAMD (for classical MD and steered/SMD runs), AMBER (via `tleap`),
Packmol (solvent/box construction), and Antechamber (small-molecule force-field
parameterisation) — plus a QM engine path-validation check (confirming XTB/ORCA/DFTB+
scratch directories and executables referenced by `$PDYNAMO3_*COMMAND` environment
variables actually exist) surfaced as a clear diagnostic dialog on project load rather
than a delayed runtime failure. DYFF, pDynamo3's own generic force field, is assignable
to an already-loaded system directly from the treeview; a structural pattern-matching
correction was added for a confirmed DYFF mistyping of guanidinium/guanidino nitrogens
(relevant to arginine side chains, among other structures) that only manifests when the
canonical ionic pattern's formal-charge precondition is not met, reusing the same
Tripos-atom-type override hook pDynamo3's own MOL2 importer uses for equivalent
corrections.

#### 3.6 Session, System, and Project Management

A Process Manager tracks background jobs with per-system filtering and multi-selection
bulk actions. Systems and visual objects can be renamed, cloned, merged, pruned, and
deleted through the treeview, with a confirmation dialog guarding deletion and dedicated
logic ensuring the (frequent) case of deleting a system while it is the sole remaining
one leaves the application's bookkeeping in a valid, well-defined "no systems loaded"
state rather than a state that later operations silently choke on. `.easy` project files
(the application's own save format) gained a saved-camera-size correction (loading a
project no longer replays the file's own stale window-size-dependent projection matrix
in place of the current session's real window size) and a compatibility shim allowing
project files written under a newer NumPy release to load correctly under this
deployment's pinned older NumPy. A "Recent Files" menu, working-folder shortcuts on the
treeview, and a "Finished!" job-completion dialog offering to pre-fill the data-import
dialog from the job's own log file round out the everyday project-management surface.

### 4. Validation Methodology

Given the interactive, GUI-driven nature of most of this software, correctness was
established primarily through direct, scripted exercising of the real application
stack — constructing a real `MainWindow`/`EasyHybridSession` instance, driving the exact
functions and, where practical, the exact GTK signal handlers a live user interaction
would invoke, and asserting on the resulting in-memory state (atom counts, bond orders,
formal charges, pDynamo `System` contents, treeview row counts) — rather than relying
solely on visual inspection or unit tests of isolated functions. For features with a
geometric or numerical correctness claim (VSEPR hydrogen placement, fragment-attachment
orientation, surface-area computation), the produced coordinates were checked
numerically against the expected geometric relationship, and in at least one case (the
fragment-orientation sign-flip bug described in Section 7) a deliberately re-introduced
"broken" code path was run side by side with the fix, under an isolated module alias, to
confirm the test itself actually discriminated correct from incorrect behaviour rather
than passing coincidentally. Every `.glade` interface-definition file touched during
development is additionally validated with a strict XML parser (`lxml`), since GTK's own
own `Gtk.Builder` and Python's standard-library `xml.etree.ElementTree` were both found
to silently accept certain malformed XML comments that the Glade Designer application
itself rejects; a full cross-check that every signal handler referenced in a `.glade`
file resolves to a real method, and every `get_object()` call resolves to a real widget
ID, was run after any interface change. Full, clean application-launch smoke tests
(construction of the main window and, where relevant, opening and closing the affected
tool windows) were run after each change to catch import-time or construction-time
regressions.

### 5. Discussion and Known Limitations

Several design decisions and open issues are worth recording explicitly. The Builder's
automatic hydrogenation and clean-up geometry are deliberately *not* a force field —
they are a cheap, purely geometric approximation intended to produce a reasonable
starting structure, not a physically accurate one; genuine energy minimisation is
available as a separate, explicit action (the DYFF "Optimize" button) precisely so the
two concerns — "give me a sensible starting guess" versus "actually relax this
structure's energy" — are not conflated. The reaction-path and conjugate-peak-refinement
tools inherit a real limitation from the underlying pDynamo3 optimiser: neither performs
an internal convergence check, so a returned path's endpoint is not, by itself, a
guarantee of having located a stationary point. A performance concern was identified,
but not yet resolved, in the bond-order perceiver used to bootstrap an "Edit in Builder"
session from an already-loaded, non-Builder-authored structure: the perception step can
be slow on a real, full-size protein, though it completes quickly on small test
molecules. Fragment attachment currently supports only a single explicit bond to the
rest of the structure (multi-point attachment is out of scope), and several
chemically reasonable but hypervalent or formally charged quick-pick fragments
(carboxylate, sulfonyl, phosphonate, nitro) are deliberately not offered yet, since the
fragment loader's present valence table assumes ordinary, non-hypervalent, neutral
main-group chemistry; extending that table is a natural, scoped next step rather than a
structural limitation of the fragment mechanism itself.

### 6. Conclusion

EasyHybrid3 combines a general-purpose QM/MM simulation library with a substantial,
independently engineered graphical and interactive layer: a full molecule builder
supporting both from-scratch construction and non-destructive editing of existing,
multi-object systems; docking, enhanced-sampling, and reaction-path workflows wired
directly into the same session model; and a broad set of structural-analysis and
visualisation tools built on a custom OpenGL renderer. The implementation history
recorded here — including the specific defects found and corrected along the way —
reflects a development process driven by continual, concrete use of the software
against real molecular systems, which is offered as part of the record precisely
because several of the more subtle fixes (index-shift hazards during cascading atom
deletion, a sign error in a 3D orientation calculation, state left inconsistent by an
earlier well-intentioned but overly destructive feature) would not have been apparent
from a purely abstract reading of the code.

---

## Part II — Detailed Changelog

Grouped by subsystem. Each entry is a short description of a distinct feature or fix;
entries are not necessarily in strict chronological order within a group.

### Molecular Builder — core editing

- Two independent Builder implementations were found in the codebase; the dead one was
  left untouched and all work targeted the live implementation.
- Fixed: delete-atom had no undo snapshot; `undo()` leaked `atom_dic_id` entries on
  every call.
- Fixed: a stale press-flag caused the render pass immediately after placing an atom to
  misinterpret it as a fresh press, starting a phantom bond-drag.
- Fixed: dragging to create a bonded atom silently discarded the new atom whenever the
  drag passed near one of the origin atom's own pre-existing hydrogens.
- Reworked: a blank Builder session no longer commits a pDynamo system until a real
  structural edit happens (previously, opening and closing the Builder without drawing
  anything still permanently created an empty system and treeview row); fixed a
  `system_liststore` duplication bug in the same pass.
- Replaced "lines" representation with "sticks" (bond-order-aware, correctly renders
  double/triple bonds as parallel cylinders) at every Builder call site.
- Added Single/Double/Triple/Aromatic bond-order selector; aromaticity tracked as a
  separate flag mirroring pDynamo3's own bond model.
- Rewrote "Clean Up" as a two-phase geometric relaxation (heavy-atom skeleton, then
  hybridisation-aware VSEPR hydrogen placement) replacing an earlier ad hoc fan
  heuristic; fixed a duplicate-hydrogen bug in the two-substituent case.
- Changed existing-bond-order editing from a Ctrl+click cycle to applying the sidebar's
  explicit selection; added, and later removed in favour of a dedicated click tool, a
  pk1/pk2 "Set Bond Order" button.
- Paired the "sticks" representation with "stick_spheres" (ball-and-stick); found and
  fixed a vismol bug rendering invisible zero-radius spheres unless a specific
  activation flag was set.
- Added drag-to-bond-by-distance (short/medium/long drag maps to triple/double/single
  bond order when the sidebar is at its default selection) with a live, colour-coded
  preview ring, and a snap-target highlight showing which existing atom a drag is about
  to bond onto.
- Added a numeric ("1"/"2"/"3") label to the bond-order preview ring, via a new minimal
  text-rendering routine (the codebase's existing hover-text routine was both dead code
  and incompatible with the current font rendering pipeline).
- Added a dedicated "Bond Order" click tool (click two atoms in sequence to change the
  order of the bond between them), superseding the pk1/pk2 button.
- Fixed a regression where a general "deactivate any representation not just rebuilt"
  sweep was silently disabling `stick_spheres` on every Builder mutation.

### Molecular Builder — fragment and structure libraries

- Added a fragment library (`.mol2` functional groups, exactly one deliberately
  under-valent atom auto-detected as the attachment point) and click-to-attach-at-a-
  hydrogen interaction; wrote a self-contained MOL2 parser after finding the codebase's
  existing one was unreachable dead code.
- Added a structure library (complete standalone molecules placed by a click, no
  hydrogen reference needed); fixed an element-resolution collision between AMBER
  atom-type codes and element symbols (`NA`/`CA`/`OS`/`NB` vs. Na/Ca/Os/Nb) found against
  a real exported `.mol2` file.
- Fixed: `.mol2` aromatic ("ar") bonds were rendered as all-single on import into the
  Structure Library; now re-Kekulised via the project's own bond-order perceiver
  whenever any bond in the file is aromatic-coded.
- Added a live, independently rendered 3D preview to both library windows; fixed sticks
  being rendered far too thin (and multi-bond separation therefore invisible) due to the
  preview's own disposable rendering session using a much smaller default stick radius
  than the main application.
- Bulk-populated the structure library with 253 real medicinal-chemistry-relevant
  molecules retrieved via the PubChem PUG REST API and processed with Open Babel, each
  validated against the real structure loader.
- Authored seven additional attachable fragments (amide, isocyanate, methoxy, hexyl,
  pentyl, phenyl, 2-furyl), each hand-verified against the fragment loader's own
  standard-valence table and confirmed to auto-detect the intended attachment atom with
  zero loader errors; four chemically requested but hypervalent/charged fragments
  (carboxylate, sulfonyl, phosphonate, nitro) were deliberately deferred, since they are
  incompatible with the loader's present non-hypervalent, neutral valence assumptions.
- Reorganised the sidebar to add eight quick-pick fragment buttons directly below the
  element selector, sharing one radio-button group with the element selector itself (so
  picking an element or a fragment are mutually exclusive, single-choice actions), with
  a "Other..." entry opening the full library window; the previously separate "Fragment"
  tool radio button was removed, its behaviour absorbed into the general "Add" tool.

### Molecular Builder — sidebar and tool selection UX

- Reorganised the sidebar into fixed-column grids to hold its width to 300 px
  regardless of the number of buttons in any one selector group; removed the pk1/pk2
  "Set Bond Order" button; moved the dihedral-angle slider out of the sidebar into its
  own window; added seven additional quick-pick elements (B, P, S, F, Cl, Br, I).
- Converted the bond-order selector to real, dot-indicator radio buttons (matching a
  user-authored glade edit that had been mistakenly reverted in an earlier pass rather
  than completed); removed the "Aromatic" bond-order option.
- Extended the "Adjust hydrogen count" switch, previously read only by the "Clean Up"
  button, to also gate every *automatic* post-edit hydrogenation call (atom placement,
  bond dragging, element replacement, bond-order changes, fragment attachment) — with
  the switch off, non-standard valence states (for example a four-bond, positively
  charged ammonium nitrogen) can be constructed by hand without the tool silently
  reverting each addition; deliberately left the explicit Protonate/Deprotonate actions
  and the Clean Up button's own independent switch unaffected, since both are meant to
  work regardless of this general switch's state.
- Fixed: deleting a heavy atom while "Adjust hydrogen count" is on left that atom's own
  hydrogens behind as disconnected, free-floating atoms; they are now removed together
  with it (deliberately without touching any *other* atom's own hydrogen count, and
  without affecting a deliberate, direct deletion of a single hydrogen atom, which must
  not be silently "undone" by this cleanup).

### Molecular Builder — geometry correctness

- Rewrote hydrogen placement to be hybridisation-aware (VSEPR-style, shared between
  automatic hydrogenation and heavy-atom skeleton relaxation) instead of a fixed-angle
  fan; fixed a planarity bug that placed a trigonal (sp2) centre's hydrogen out of plane
  whenever its two existing substituents were not exactly 120° apart (in practice, for
  almost any real structure).
- Fixed a genuine 180° orientation error in fragment attachment: the rotation aligning a
  newly attached fragment with its bond to the rest of the molecule used the wrong sign
  for the target direction vector, causing every attached fragment to fold back toward
  the parent structure instead of extending away from it, independent of where on the
  molecule the attachment occurred. The correct sign was derived by cross-checking
  against the already-correct, structurally identical use of the same underlying
  direction-completion routine in ordinary hydrogen placement, and confirmed with a
  side-by-side comparison against a deliberately reintroduced copy of the old, incorrect
  sign.

### Molecular Builder — protonation state

- Added a `formal_charge` attribute to the visual atom model, deliberately kept
  separate from the pre-existing `charge` attribute (an MM partial charge populated only
  once a force field is assigned, an unrelated quantity already used elsewhere in the
  codebase).
- Extended the automatic-hydrogenation target-valence calculation to account for formal
  charge (a `+1` charge raises target valence by one bond, matching one fewer available
  lone pair; a `-1` charge lowers it by one, symmetrically), making protonation and
  deprotonation a matter of changing the charge by one unit and letting the existing
  hydrogenation machinery place or remove the corresponding proton.
- Added "Protonate"/"Deprotonate" actions and an editable "Charge" column to the "Atom
  Types" window, operating on the current multi-row selection; fixed a real
  correctness bug where a naive batch implementation could target the wrong atom
  partway through a multi-atom protonation/deprotonation, since removing one atom's
  hydrogen shifts the internal indices of every atom numbered above it.
- Propagated formal charge into the real pDynamo `Atom.formalCharge` field whenever the
  underlying system is (re)built, fixing a previously accepted, documented limitation
  ("no formal-charge UI exists") that had already been identified as the root cause of a
  separate DYFF force-field mistyping issue (see the guanidinium correction below).

### Molecular Builder — editing an already-loaded system

- Added the ability to edit an already-loaded system's visual object directly, via a
  disposable clone; the user selects which visual object to edit when a system has more
  than one, and is warned that an existing force field or QM model assignment will be
  lost by editing (since editing necessarily rebuilds the underlying pDynamo `System`
  from scratch).
- Implemented the fold-back-versus-new-system decision: an edit that preserves atom
  count is folded back onto the original system in place; an edit that changes atom
  count instead becomes its own new, independent system, since folding a changed atom
  count back in place would desynchronise other atom-index-keyed application state.
- Fixed: a system-editing session, if left open, could silently pick a molecular-surface
  child object as its "original" target purely by dictionary-iteration-order luck; fixed
  by excluding surface objects from the fallback lookup, matching the exclusion already
  applied elsewhere.
- Fixed: after cloning a normally loaded system, the clone could inherit other active
  representations (for example a protein's cartoon ribbon) that the Builder does not
  keep synchronised, causing an index-out-of-bounds crash on the very first click;
  those representations are now deactivated immediately on entering edit mode.
- Fixed a real, user-reported non-destructive-editing regression: an earlier
  implementation unconditionally discarded every *other* visual object belonging to the
  same system before editing began, regardless of whether the session would eventually
  fold back in place or become its own new system — corrupting the original system even
  in the "new system" case, where it is supposed to remain completely untouched. The
  discard step was removed entirely; the specific visual object being edited is now
  tracked by a direct object reference rather than an ambiguous "first object matching
  this system" lookup, which remains correct even when a system genuinely has more than
  one associated visual object at fold-back time.
- Added reversible visibility hiding: every other currently visible object (including
  the very object now being edited, since the Builder works on a separate clone of it)
  is temporarily deactivated for the duration of an editing session and restored,
  unconditionally, when the session ends — implemented as a pure visibility toggle,
  never touching any object's underlying data, and tracking exactly which objects it
  hid so an object that was already hidden for an unrelated reason is never incorrectly
  reactivated.
- Fixed a "New (blank) system" regression that only worked the first time: closing the
  Builder after a blank session that had received a real edit failed to reset the
  session's internal target-object pointer, so the next "New" request silently resumed
  the previous, already-promoted system instead of starting fresh; separately fixed
  requesting "New" while a Builder session was already open (without explicitly
  closing it first) silently continuing to edit the existing session instead of starting
  a new one.
- Fixed a deeper instance of the same class of bug: deleting the sole remaining loaded
  system leaves a documented, intentional placeholder value in the application's system
  table; a still-open Builder session whose underlying system had just been deleted
  could, on its next edit, mistake that placeholder for a live system and register a
  new system missing required one-time bookkeeping attributes, silently and
  permanently breaking the systems tree view's ability to render *any* further system
  from that point on. Both the tree view's own rendering loop and the session-
  registration logic were fixed to correctly treat the placeholder as "no system here."

### Molecular Builder — auxiliary tools

- Added a periodic-table picker (full main-group and transition-metal range, UFF-derived
  per-element data reused from the same source pDynamo3's own generic force field
  cites) alongside the sidebar's quick-pick element buttons.
- Added an "Atom Types" window showing and allowing manual correction of each selected
  atom's live-perceived DYFF force-field type, without requiring a force field to
  already be assigned.
- Fixed a systematic half-character-width centring offset in the ring/label overlay used
  by several Builder tools (bond order, dihedral, atom-type previews).
- Added a one-click "Optimize (DYFF)" action running a bounded, real force-field
  minimisation on a disposable scratch copy of the working structure.
- Added a "Transform Selection" window (six sliders: translate/rotate about X/Y/Z,
  acting on the current selection about its own centroid).
- Added a dedicated dihedral-angle rotation window (pick four atoms, then drive one side
  of the bond with a slider), extracted from an earlier inline sidebar slider.
- Suppressed the application's general-purpose shift-click/shift-drag viewing selection
  — both its rendering and its underlying state — while the Builder's own editing mode
  is active, and restored it automatically when editing is paused or the Builder is
  closed, so the two independent selection mechanisms never interfere with one another.

### Molecular Docking

- Added an AutoDock Vina docking front end (ported from an existing standalone tool);
  fixed a shared-`GtkAdjustment` bug, an Open Babel duplicate-atom-name import bug, and
  a crash loading a headerless PDB.
- Added an AutoDock-GPU docking front end, reusing the Vina tool's own code; fixed an
  AutoGrid4 crash caused by an excessively long grid-file path.
- Reworked pose import to write directly into `VisMol` objects backed by a real,
  in-memory pDynamo `System`, removing an intermediate temporary-PDBQT round trip, and
  added a per-object mechanism for pinning which trajectory frame a pose-derived object
  displays.

### Enhanced Sampling and Reaction-Coordinate Methods

- Added two-dimensional potential-energy-surface scanning with dihedral reaction-
  coordinate support, matching the pre-existing one-dimensional case.
- Merged the "simple" (fixed-shape reaction coordinate) and "advanced" (arbitrary
  weighted-distance-sum reaction coordinate) potential-energy-surface scan windows into
  one, behind a mode toggle, with the computational backend unchanged.
- Applied the same simple/advanced merge to the Umbrella Sampling window, and added
  genuine dihedral-restraint support (correct 360° periodicity) rather than
  distance-only restraints; fixed a dead/broken parameter-restoration code path found
  in the same pass.
- Added an "Energy Refinement" tool that measures (rather than scans) a reaction
  coordinate along an existing trajectory, sharing the advanced-mode weighted-sum
  reaction-coordinate logic; fixed a log-file-type recognition bug that silently
  prevented its own log output from being re-imported for plotting.
- Added automatic reaction-coordinate auto-detection from a chosen results folder's
  log file (for both Umbrella Sampling and Energy Refinement), recognising three
  historical log formats.
- Added an interactive statusbar readout and curve-identification-by-click to the WHAM
  post-processing windows, and automatic PNG export of the resulting histograms and
  free-energy profile on job completion.
- Added a conjugate peak refinement (CPR) reaction-path tool, built by generalising the
  nudged-elastic-band tool's own initial-path construction.

### Trajectory and Structural Analysis

- Implemented a from-scratch cartoon (ribbon) representation, including a substantial
  vectorised performance pass and an optional dynamic-secondary-structure caching mode.
- Implemented a general molecular-surface tool (van der Waals, solvent-accessible, and
  solvent-excluded surfaces), coloured by the viewer's own existing per-atom colours,
  including a mode colouring only carbon atoms by a system's custom reference colour;
  added terminal-only total surface-area reporting, structured to support a future
  partial-surface query.
- Added radial distribution function, RMSF, and Ramachandran analysis windows.
- Fixed three bugs in the infrared-spectrum tool's underlying dipole-moment calculation
  (`QCModelXTB.DipoleMoment`), and documented that the system-level dipole moment does
  not sum contributions across mixed MM/QC regions.

### QM/MM Setup and External Engine Integration

- Added setup front ends for NAMD, AMBER (`tleap`), Packmol, and Antechamber.
- Added a QM engine path-validation check (XTB/ORCA/DFTB+ scratch directories and
  executables), surfaced as a diagnostic dialog on project load.
- Added the ability to assign pDynamo3's DYFF generic force field to an already-loaded
  system from the treeview; fixed a parameter-set auto-discovery bug and a confirmed
  DYFF mistyping of guanidinium/guanidino nitrogens (relevant to arginine side chains),
  corrected via the same Tripos-atom-type override hook pDynamo3's own MOL2 importer
  uses.

### Session, System, and Project Management

- Added a Process Manager with per-system filtering and multi-selection bulk actions.
- Added system/visual-object renaming, and a confirmation dialog guarding system
  deletion.
- Fixed a `.easy` project file loading bug where the saved camera projection used the
  file's own stale window-size metadata instead of the current session's real window
  size, and added a NumPy pickle-compatibility shim for project files written under a
  newer NumPy release than this deployment's pinned version.
- Added a "Recent Files" menu, working-folder shortcuts on the systems tree view, and a
  job-completion dialog offering to pre-fill the data-import dialog from the
  completed job's own log file.

---

*This document was generated as an internal implementation report accompanying ongoing
development of EasyHybrid3. It reflects the state of the codebase at the time of
writing and should be revised as the software continues to evolve.*
