# Metalloprotein test set

Well-studied metalloproteins from the RCSB PDB (https://www.rcsb.org), for testing the
DYFF preparation without OPLS (`src/pdynamo/dyff_prep.py`, Extras > Prepare DYFF System...).
Download / refresh: `python3 fetch_metalloproteins.py` (writes `manifest.json` with title,
method, resolution and the metal sites read from the LINK records or from distances).

| PDB | Protein | Method | Resolution (A) | Site (literature) | Why it is in the set |
|---|---|---|---|---|---|
| 1ZNF | Zinc finger Xfin-31 (NMR) | SOLUTION NMR | - | Zn Cys2His2 (structural) | classic structural Zn site; NMR model, caps ACE/NH2 |
| 3KS3 | Carbonic anhydrase II, human (0.9 A) | X-RAY DIFFRA | 0.9 | Zn His3 + OH-/H2O (catalytic) | the reference catalytic Zn site |
| 8TLN | Thermolysin | X-RAY DIFFRA | 1.6 | Zn His2Glu + H2O; four Ca sites | Zn protease with structural Ca |
| 4INS | Insulin hexamer (2Zn) | X-RAY DIFFRA | 1.5 | Zn His3 x2 (on the 3-fold axis) | Zn at symmetry axis, multi-chain |
| 2C9V | Cu,Zn superoxide dismutase 1, human | X-RAY DIFFRA | 1.07 | Cu His4 (bridging His); Zn His3Asp | binuclear Cu/Zn with bridging imidazolate |
| 4AZU | Azurin, Pseudomonas aeruginosa | X-RAY DIFFRA | 1.9 | Cu(II) type 1: His2 Cys Met (+Gly O) | blue copper, open-shell Cu(II) |
| 1PLC | Plastocyanin, poplar (1.33 A) | X-RAY DIFFRA | 1.33 | Cu type 1: His2 Cys Met | blue copper reference |
| 1IRO | Rubredoxin, Clostridium pasteurianum (1.1 A) | X-RAY DIFFRA | 1.1 | Fe Cys4 (high-spin Fe(II)/Fe(III)) | simplest Fe-S site, open shell |
| 1FDN | Ferredoxin 2[4Fe-4S], C. acidurici | X-RAY DIFFRA | 1.84 | two [4Fe-4S] clusters, Cys4 each | iron-sulfur clusters (cluster residues) |
| 1A6M | Myoglobin, sperm whale, oxy (1.0 A) | X-RAY DIFFRA | 1.0 | heme Fe + proximal His + O2 | heme cofactor (porphyrin as one residue) |
| 1HRC | Cytochrome c, horse heart | X-RAY DIFFRA | 1.9 | heme c, Fe His/Met, thioether links to Cys | covalently attached heme |
| 1CLL | Calmodulin, human (Ca4) | X-RAY DIFFRA | 1.7 | four EF-hand Ca sites (Asp/Glu carboxylates, O) | hard ion, O-only coordination |
| 5P21 | H-Ras p21 with GppNHp and Mg | X-RAY DIFFRA | 1.35 | Mg Ser/Thr + phosphate oxygens + waters | Mg with nucleotide |
| 1UBQ | Ubiquitin, human (no metal) | X-RAY DIFFRA | 1.8 | negative control: no metal | checks that nothing is detected |

Notes
- 1ZNF is an NMR ensemble (37 models; the first is used) with ACE/NH2 caps bonded to the
  chain: the preparation detects them as covalently bound and removes them.
- 4INS: each Zn is coordinated by His B10 of three chains related by crystal symmetry; the
  asymmetric unit holds one His per Zn, so the site is incomplete (treated as ionic, with a
  warning) unless the biological assembly is built.
- 1FDN, 1A6M, 1HRC: metal-containing cofactors (Fe-S clusters, heme) -- not supported yet by
  the DYFF preparation (reported).
- 1UBQ: negative control, no metal.
