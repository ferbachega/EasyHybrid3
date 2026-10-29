# Tutorial input files

The input files of the tutorials and case studies of the EasyHybrid site (`docs/site/`), one folder
per tutorial, each with the tutorial itself as a PDF. Paths in the tutorials are relative to the
EasyHybrid folder. The PDFs are written by `docs/site/build_site.py` (front matter key `examples`).

## Tutorials

| # | Tutorial (`docs/site/`) | Files |
|---|---|---|
| 1 | Blocked alanine: optimization and a first scan (`bala_tutorial.md`) | `01_bala/bALA_not_opt.mol` |
| 2 | A simple SN2 reaction with a QC method (`sn2_tutorial.md`) | `02_sn2/SN2_tutorial.xyz` |
| 3 | Diels–Alder reaction: a 2D energy surface (`diels_alder_tutorial.md`) | `03_diels_alder/diels_alder.mol2` |
| 4 | Molecular dynamics basics (`md_basics_tutorial.md`) | `04_md_basics/1UAO_OPLS_folded_solvated.pkl` |
| 5 | Conformational change with NEB: cyclohexane (`neb_cyclohexane_tutorial.md`) | `05_neb_cyclohexane/cyclohexane_chair.mol`, `cyclohexane_twistboat.xyz` |
| 6 | A first QC/MM system (`qcmm_tutorial.md`) | `06_qcmm/1UAO_OPLS.pkl` |
| 7 | Umbrella sampling of blocked alanine (`bala_umbrella_tutorial.md`) | `07_bala_umbrella/` (PDF only): the scan `bALA_scan_OH.ptGeo` made in tutorial 1 |
| 8 | QC/MM free energy of the SN2 step in LinB (`linb_sn2_tutorial.md`) | `08_linb_sn2/LinB.yaml` |
| 9 | Metalloprotein with DYFF (`pages/metalloprotein_tutorial.html`) | `09_metalloprotein/1ZNF.pdb` |
| 10 | Radial distribution functions (`rdf_tutorial.md`) | `10_rdf/` (PDF only); the system is `04_md_basics/1UAO_OPLS_folded_solvated.pkl` |
| 11 | Checking and exploring umbrella sampling runs (`umbrella_browser_tutorial.md`) | `11_umbrella_browser_sn2/`: system, scans and finished US runs (1D, 1D stiff, 2D) |
| 12 | EasyHybrid with Claude Code and a local model (`ai_assistants_tutorial.md`) | `12_ai_assistants/` (PDF only); the system is `01_bala/bALA_not_opt.mol` |

## Case studies

| # | Case study | Files |
|---|---|---|
| 1 | Triosephosphate isomerase (`tim_case_study.md`) | `case_studies/1_tim/TIM_DHAP_AMBER_AM1.pkl` |
| 2 | Chorismate mutase (1): scan (`cm_scan_case_study.md`) | `case_studies/2-4_chorismate_mutase/chorismate_mutase_AMBER_PM6.pkl` |
| 3 | Chorismate mutase (2): NEB (`cm_neb_case_study.md`) | the scan `CM_scan.ptGeo` made in case study 2 |
| 4 | Chorismate mutase (3): umbrella sampling (`cm_us_case_study.md`) | the scan `CM_scan.ptGeo` made in case study 2 |
| 5 | β-Phosphoglucomutase (`pgmase_case_study.md`) | `case_studies/5_pgmase/phosphoglucomutase_AMBER_AM1dPhoT.pkl` |

## User guide

- Transition state search (`ts_search_manual.md`): `05_neb_cyclohexane/cyclohexane_chair.mol` (the
  same file as `$PDYNAMO3_HOME/examples/book/data/mol/cyclohexane_chair.mol`).

Note: `*.pkl` is in `.gitignore`; add the `.pkl` files here with `git add -f`.
