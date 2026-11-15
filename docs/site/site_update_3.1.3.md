# EasyHybrid Google Site: content for the 3.1.3 update

Draft of 2026-10-16. Each section below is one page, or one block of a page, on
https://sites.google.com/view/easyhybrid. Lines marked **[CONFIRM]** need a decision before
publishing.

Out of scope for now, by request: the Builder, the NAMD interface, AutoDock Vina, AutoDock-GPU
and the other Extras tools.

---

## 1. Home: new block "What's new in 3.1.3" (under the feature list)

**What's new in EasyHybrid 3.1.2 – 3.1.3** (September–October 2026)

- **IR spectra from normal modes.** Infrared intensities, a symmetry label for each
  vibrational mode, and a broadened IR spectrum plot.
- **Molecular surfaces.** Van der Waals, solvent-accessible and solvent-excluded (Connolly)
  surfaces for any molecule, colored by atom, chain or carbon color.
- **New cartoon representation.** Written from scratch, with its own secondary-structure
  assignment, several color schemes and transparency. It follows trajectories.
- **Interface themes and depth cueing.** Five color themes for the interface, plus a
  configurable fog in the 3D view.
- Fixes for macOS graphics, and many smaller fixes.

[See the release notes and download →](Downloads)

---

## 2. Downloads (replaces the current page)

### Current version: 3.1.3 (October 5, 2026)

**[CONFIRM]** On GitHub, `main` is still at 3.1.2; 3.1.3 is only on `dev`. Either merge `dev`
into `main` before publishing this page, or present 3.1.3 as the development version.

#### Requirements

- **pDynamo3.** Required. A partner project distributed separately:
  https://github.com/pdynamo/pDynamo3
- **Python 3** with: numpy, PyOpenGL, Cython, freetype-py, pycairo, PyGObject, Pillow, psutil.
  All of them are listed in `requirements.txt`:

      pip install -r requirements.txt --user

- **GTK 3** typelib. This is an operating-system package, not a Python one. On Debian/Ubuntu:

      sudo apt-get install gir1.2-gtk-3.0

- **macOS only:** the GLFW library, for example `mamba install -c conda-forge glfw`.
- **Optional:** xTB (GFN1/GFN2-xTB), ORCA and other external QC programs, set up in
  Preferences.

`install.py` checks all of the above and reports anything that is missing.

#### Platforms

Linux, macOS and Windows (experimental).

#### Stable release

    git clone --recurse-submodules https://github.com/ferbachega/EasyHybrid3.git
    cd EasyHybrid3
    python3 ./install.py
    ./easyhybrid.py

Check that `src/graphics_engine` exists after cloning. If it is missing, the submodule was not
downloaded: run `git submodule update --init`.

#### Development builds

Beta versions, updated almost weekly. Use at your own risk.

    git clone --recurse-submodules -b dev https://github.com/ferbachega/EasyHybrid3.git

#### Release notes: 3.1.2 and 3.1.3

(IR intensities and spectra arrived in 3.1.2; everything else in this list is from 3.1.3.)

**Analysis**
- IR intensities in *Simulate › Normal Modes* (option *Compute IR intensities*).
- Symmetry (irreducible representation) of each vibrational mode, and an *IR Intensity*
  column in *Analysis › Normal Modes*.
- *Plot IR Spectrum* button: stick spectrum plus a Lorentzian-broadened curve.

**Visualization**
- Molecular surfaces: VWS, SAS and SES, with probe radius, grid spacing, Taubin smoothing,
  colors by atom, chain or carbon, and the surface area. Each surface is its own object in
  the tree. Surfaces can be generated for every frame of a trajectory.
- New cartoon representation:
  - its own secondary-structure assignment (Cα-based, P-SEA-like), optionally reassigned on
    every frame;
  - colors by secondary structure, rainbow, chain, atom color, or a single color;
  - transparency, and click-to-select on the cartoon;
  - adjustable helix, strand, arrow and coil dimensions (*Preferences › Viewer › Representations › Cartoon*).
- Fog (depth cue): on/off, intensity, automatic or custom range (*Preferences › Viewer ›
  General › Fog*).

**Interface**
- Interface themes (*Preferences › Interface*, *Interface theme*), applied live: Scientific
  Ribbon, Graphite, Nord, Slate (light) and High contrast, plus the system GTK theme.
- *File › Open Recent*.
- macOS rendering and trackpad fixes.

#### License

EasyHybrid and its graphics engine are distributed under the GNU GPLv3. The core and graphics
libraries are under continuous development.

#### Acknowledgment

(keep the current paragraph thanking Dr. Ahmad Shabir)

---

## 3. Citations (replaces the current page)

If EasyHybrid contributes to your work, please cite:

**EasyHybrid**

> Bachega, J. F. R.; Hagen, G.; Sequeiros-Borja, C.; Nikklas, K.; Chahine, J.;
> Timmers, L. F. M. S.; Field, M. J. EasyHybrid: An Interactive Graphical Environment for
> Quantum, Classical and Hybrid Simulations with pDynamo3. *J. Chem. Inf. Model.* **2026**,
> *66* (3), 1286–1292. https://doi.org/10.1021/acs.jcim.5c02047

(Published online on November 11, 2025. Selected for the cover of the *Journal of Chemical
Information and Modeling*. Keep the current cover image.)

**pDynamo**

> Field, M. J. pDynamo3: A Modern Framework for Molecular Simulations Using Hybrid Quantum
> and Classical Methods. *J. Chem. Inf. Model.* **2022**, *62* (23), 5849–5854.
> https://doi.org/10.1021/acs.jcim.2c01239

> Field, M. J. The pDynamo Library for Molecular Simulations Using Hybrid Quantum Mechanical
> and Molecular Mechanical Potentials. *J. Chem. Theory Comput.* **2008**, *4* (7),
> 1151–1161. https://doi.org/10.1021/ct800092p

> Field, M. J. *A Practical Introduction to the Simulation of Molecular Systems*, 2nd ed.;
> Cambridge University Press, 2007.

**BibTeX**

```bibtex
@article{Bachega2026EasyHybrid,
  author  = {Bachega, Jos{\'e} Fernando R. and Hagen, Gabriel and Sequeiros-Borja, Carlos and
             Nikklas, Kevin and Chahine, Jorge and Timmers, Luis Fernando S. M. and Field, Martin J.},
  title   = {{EasyHybrid}: An Interactive Graphical Environment for Quantum, Classical and
             Hybrid Simulations with {pDynamo3}},
  journal = {Journal of Chemical Information and Modeling},
  year    = {2026},
  volume  = {66},
  number  = {3},
  pages   = {1286--1292},
  doi     = {10.1021/acs.jcim.5c02047}
}

@article{Field2022pDynamo3,
  author  = {Field, Martin J.},
  title   = {{pDynamo3}: A Modern Framework for Molecular Simulations Using Hybrid Quantum
             and Classical Methods},
  journal = {Journal of Chemical Information and Modeling},
  year    = {2022},
  volume  = {62},
  number  = {23},
  pages   = {5849--5854},
  doi     = {10.1021/acs.jcim.2c01239}
}

@article{Field2008pDynamo,
  author  = {Field, Martin J.},
  title   = {The {pDynamo} Library for Molecular Simulations Using Hybrid Quantum Mechanical
             and Molecular Mechanical Potentials},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2008},
  volume  = {4},
  number  = {7},
  pages   = {1151--1161},
  doi     = {10.1021/ct800092p}
}
```

**[CONFIRM]** The co-authors' first names in the BibTeX (Gabriel, Carlos, Kevin, Jorge,
Luis Fernando) were not checked against the paper. Use initials if unsure.

---

## 4. User Guide: new pages

Links on the User Guide index:
- *Interface basics*: add **Interface themes and fog**.
- *Interface basics*: add **Cartoon representation**, under *Colors and representations*.
- *Analysis*: link the existing **Surfaces** item to **Molecular surfaces**.
- *Analysis*: link the existing **Normal modes** item to **Normal modes and IR spectra**.

### 4.1 Normal modes and IR spectra

EasyHybrid computes vibrational normal modes for any system with an energy model: QC, MM or
hybrid QC/MM. Since version 3.1.2 it can also compute infrared intensities, label each mode by
its symmetry, and plot the IR spectrum.

**1. Run the normal-mode calculation**

1. Optimize the geometry first. Frequencies are only meaningful at a minimum; imaginary
   frequencies mean the structure is not one.
2. Open *Simulate › Normal Modes*.
3. Tick **Compute IR intensities**. It is off by default because it is expensive:
   - the dipole moment is differentiated numerically;
   - this needs about six extra energy-and-dipole evaluations per free atom.
4. Click **Run**.

**2. Inspect the modes**

Open *Analysis › Normal Modes*, choose the system and import the result. The table lists, for
each mode:
- the frequency (cm⁻¹);
- the **Symmetry**: the irreducible representation of the molecule's point group;
- the **IR Intensity** (km/mol). It shows `?` when the run did not compute intensities.

Tick a mode to animate it in the 3D view.

**3. Plot the spectrum**

Click **Plot IR Spectrum**. The plot shows:
- the stick spectrum (grey);
- a curve broadened with Lorentzian lines of 20 cm⁻¹ FWHM (orange), as an experimental
  spectrum would look.

Imaginary modes, and modes below 50 cm⁻¹ (residual translations and rotations), are left out.

**Notes and limits**

- **Hybrid QC/MM:** pDynamo's dipole moment only includes the QC region. The intensities then
  describe the QC region only, and the log prints a warning.
- **MM models need charges.** A force field without electrostatics gives zero intensities;
  for DYFF or UFF, assign charges first (for example CM5).
- **QC engines:** semi-empirical (MNDO family) and ORCA work with pDynamo3 as distributed.
  **[CONFIRM]** xTB needs the updated xTB interface in pDynamo3 (dipole read from
  `xtbout.json`); say where users get it, or leave xTB out of the list.
- Symmetry labels come from pDynamo's point-group perception. A slightly distorted geometry
  can lower the detected symmetry.

### 4.2 Molecular surfaces

EasyHybrid draws molecular surfaces for any object, protein or small molecule, independent of
the QC surfaces (orbitals, density, electrostatic potential).

**Surface types**

| Type | What it is |
|---|---|
| **Van der Waals (VWS)** | The union of the atomic van der Waals spheres. |
| **Solvent-Accessible (SAS)** | Traced by the center of a probe sphere rolling over the molecule: the spheres grown by the probe radius. |
| **Solvent-Excluded / Connolly (SES)** | Traced by the surface of the probe: the region the solvent cannot reach. The usual "molecular surface". |

**Creating a surface**

1. Right-click the object in the tree, then choose *Representation › New Surface…*.
2. In **Surface Setup**, choose:
   - **Surface type**: VWS, SAS or SES (default SES).
   - **Probe radius** (SAS/SES only): 1.4 Å for water.
   - **Grid spacing**: default 0.8 Å. Smaller values give a finer surface but are slower.
   - **Color by**: *Atom colors*, *Chain*, or *Carbon color (system)*.
   - **Smooth (Taubin)**, optional, with the number of iterations.
3. Click **Generate**.

**How surfaces behave**

- Each surface is a new object, nested under its molecule in the tree. You can have several
  at once, even of the same type: SAS and SES, or different probe radii.
- For a trajectory, a surface is generated for every frame (in parallel), and it follows the
  frame slider.
- The surface area (Å²) is printed in the terminal when the surface is generated.

**Changing a surface afterwards**

Right-click the surface object, then choose *Setup*. Options:
- **Wireframe**;
- **Opacity**;
- **Smooth shading**;
- **Decimate**, which merges vertices closer than a given distance to make the mesh lighter;
- **Smooth (Taubin)**.

Default values for new surfaces are in *Edit › Preferences › Viewer › Representations › Surfaces*.

### 4.3 Cartoon representation

The cartoon was rewritten in 3.1.3. It is EasyHybrid's own implementation:
- helices and strands are flat ribbons, with arrowheads on strands;
- coils are round tubes;
- secondary structure is assigned from the Cα trace (a P-SEA-like method), so no DSSP is
  needed.

**Showing the cartoon**

- One object: right-click it in the tree, then *Representation › Cartoon*.
- A whole system: right-click it, then *Representation › Cartoon*.
- Breaks between chains and gaps in the chain are respected.

**Colors: *Representation › Cartoon Setup…***

- *Secondary structure* (default; helix, strand and coil colors are set in Preferences).
- *Rainbow* (N → C, per chain).
- *Chain*.
- *Atom color (Cα)*.
- *Single color*.

**Fine tuning: *Edit › Preferences › Viewer › Representations › Cartoon***

- **Ribbon thickness:**
  - helix width and height, strand width and height;
  - arrowhead width and taper (in residues);
  - coil radius;
  - ribbon detail (points per residue);
  - strand smoothing.
- **Secondary structure assignment:**
  - helix and strand tolerances;
  - minimum helix and strand lengths;
  - **Dynamic**: reassigns helices and strands on every trajectory frame. It is off by
    default, which is faster; the cartoon still follows the coordinates.
- **Colors:**
  - helix, strand and coil colors;
  - **Transparency (%)**.

Clicking the cartoon selects the residue under the pointer.

### 4.4 Interface themes and fog

**Interface themes**

*Edit › Preferences › Interface* (Startup section, *Interface theme*). The theme is applied immediately; click
*Apply and Save* to keep it.

| Theme | |
|---|---|
| Default | the system GTK theme |
| Scientific Ribbon | dark navy |
| Graphite | dark, orange accent |
| Nord | dark, frost accent |
| Slate | light, blue accent |
| High contrast | black / yellow |

Dark themes also lighten the toolbar icons so they stay visible.

**Fog (depth cue)**

*Edit › Preferences › Viewer › General › Fog*. Fog fades distant atoms into the background and
improves the sense of depth.

- **Enable fog**: turn it on or off.
- **Intensity (%)**: 100% is the classic look; lower values give a subtler fog.
- **Automatic range**: the fog follows zoom and clipping. It starts either at the far plane
  or a few Å before the rotation center (*Offset before center*).
- **Custom range**: start and end set as a percentage of the visible depth. They still follow
  zoom.

Changes are applied immediately.

---

## 5. DYFF manual: hidden page

Publish as a page that is **hidden from navigation** (Google Sites: page ⋮ › *Hide from
navigation*). Note that a hidden page is still public for anyone who has its URL.

**[CONFIRM]** The manual describes the current development code, not 3.1.3: the DYFF residue
library, xTB refinement and Prepare DYFF System came after October 5. That is a further reason
to keep it hidden until the next release.

The full text is in `docs/site/dyff_manual.md`.
