# Tinker OPLS parameter files

Copied from the user's reference files (2026-10-04):

| File | Force field | Notes |
|---|---|---|
| `oplsaal.prm` | OPLS-AA/L (Kaminski, Friesner, Tirado-Rives, Jorgensen, J. Phys. Chem. B 2001, 105, 6474) | vdW indexed by CLASS; protein torsions reparametrized |
| `oplsaa08.prm` | OPLS-AA 2008 (BOSS 4.8, W. L. Jorgensen) | vdW indexed by TYPE; 906 atom types (organic molecules, ions, nucleic acids) |

Read by `pdynamo/opls/tinker_prm.py`; merged into the class-level database of
`pdynamo/opls/parameters.py` with precedence pDynamo `opls/protein` > `oplsaal` > `oplsaa08`.
Coverage and every conflict between the three sources: `../opls_database_report.md`
(regenerate with `parameters.database_report_markdown()`).

Before distributing EasyHybrid with these files, check the redistribution terms of the
Tinker parameter files.
