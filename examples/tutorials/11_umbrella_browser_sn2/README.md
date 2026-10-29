# Umbrella sampling example: SN2 Cl⁻ + CH₃Br (AM1)

Data of the tutorial *Checking and exploring umbrella sampling runs*
(`docs/site/umbrella_browser_tutorial.md`, Analysis → Umbrella Sampling Browser).

| Path | Content |
|---|---|
| `sn2_cl_ch3br.xyz` | starting structure (reactant complex); AM1, charge −1, all atoms QC |
| `SN2_scan_1D.ptGeo/` | relaxed scan of d(C–Br) − d(C–Cl), −1.4 to +1.6 Å, 31 frames |
| `SN2_scan_2D.ptGeo/` | relaxed 2D scan, d(C–Br) 1.8–3.2 Å × d(C–Cl) 1.7–3.1 Å, 15 × 15 frames |
| `SN2_US_1D/` | 1D umbrella sampling, k = 600 kJ mol⁻¹ Å⁻², 2 ps equilibration + 10 ps per window, WHAM in `wham/` |
| `SN2_US_1D_stiff/` | the same with k = 4000: windows that do not overlap |
| `SN2_US_2D/` | 2D umbrella sampling, k = 600, 1 ps + 4 ps per window, WHAM in `wham/` |

Langevin dynamics, 300 K, 1 fs, collision frequency 25 ps⁻¹. Atom indices (0-based): C 0,
H 1–3, Br 4, Cl 5.
