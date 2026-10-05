"""
Benchmark: OPLS classes and CM5 charges of the simplified ligand parametrization
(src/pdynamo/opls/ligand.py) against the OPLS-AA 2008 reference types.

14 small molecules with explicit OPLS-AA 2008 types are built from SMILES (Open
Babel --gen3d), parametrized with GFN1-xTB CM5 charges, and compared atom by atom
with the reference type (class) and charge. Reports the class agreement, the
charge MAE/RMSD for several CM5 scale factors and the least-squares optimum.

Needs: pDynamo3 (PDYNAMO3_HOME or ~/programs/pDynamo3), obabel, xtb.
Usage: python3 tests/benchmarks/opls_ligand_charges.py [--json out.json]
Result of 2026-10-04: classes 133/134, MAE 0.054 (s=1.00), 0.082 (1.20), optimum s = 0.835.
"""
import json, os, subprocess, sys, tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path += [os.environ.get("PDYNAMO3_HOME", os.path.expanduser("~/programs/pDynamo3")),
             os.path.join(REPO, "src"), os.path.join(REPO, "src", "graphics_engine", "src")]

from pdynamo.opls import ligand as L                    # noqa: E402
from pdynamo.opls.parameters import TINKER_FILES        # noqa: E402
from pdynamo.opls.tinker_prm import TinkerPRM           # noqa: E402

# name: ( SMILES, reference types of the heavy atoms in SMILES order,
#         { heavy atom index: type of its hydrogens }, default hydrogen type )
REFERENCE = {
    "methanol":          ("CO",             [99, 96],                  {0: 98, 1: 97}, 85),
    "ethanol":           ("CCO",            [80, 99, 96],              {2: 97}, 85),
    "dimethyl ether":    ("COC",            [123, 122, 123],           {0: 127, 2: 127}, 85),
    "acetone":           ("CC(=O)C",        [80, 222, 223, 80],        {}, 85),
    "acetic acid":       ("CC(=O)O",        [80, 209, 210, 211],       {3: 212}, 85),
    "N-methylacetamide": ("CC(=O)NC",       [80, 177, 178, 180, 184],  {3: 183}, 85),
    "methylamine":       ("CN",             [733, 730],                {1: 739}, 85),
    "ethanethiol":       ("CCS",            [80, 148, 142],            {2: 146}, 85),
    "nitromethane":      ("C[N+](=O)[O-]",  [703, 701, 702, 702],      {0: 704}, 85),
    "benzene":           ("c1ccccc1",       [90] * 6,                  {}, 91),
    "phenol":            ("Oc1ccccc1",      [109, 108] + [90] * 5,     {0: 110}, 91),
    "chlorobenzene":     ("Clc1ccccc1",     [206, 205] + [90] * 5,     {}, 91),
    "imidazole":         ("[nH]1cncc1",     [498, 499, 500, 501, 502], {0: 503, 1: 504, 3: 505, 4: 506}, 91),
    "propane":           ("CCC",            [80, 81, 80],              {}, 85),
}


def parametrize(smiles):
    scratch = tempfile.mkdtemp()
    subprocess.run(["obabel", "-:" + smiles, "--gen3d", "-h", "-O", os.path.join(scratch, "m.mol2")], capture_output=True)
    atoms, bonds = L._parse_mol2(os.path.join(scratch, "m.mol2"))
    elements = [a[1] for a in atoms]
    kekule, _ = L._kekulize(elements, bonds)
    named = [("%s%d" % (e, k + 1), L.ELEMENT_Z[e], a[2]) for k, (a, e) in enumerate(zip(atoms, elements))]
    return L.parametrize_residue("LIG", named, bonds=kekule, cm5_scale=1.0)       # unscaled CM5


def main():
    prm = TinkerPRM(dict(TINKER_FILES)["oplsaa08"])
    rows, reference_q, raw_q = [], [], []
    agree_total = 0
    for name, (smiles, heavy, h_types, h_default) in REFERENCE.items():
        lig = parametrize(smiles)
        neighbors = {}
        for i, j, _ in lig.bonds:
            neighbors.setdefault(i, []).append(j); neighbors.setdefault(j, []).append(i)
        types = [heavy[k] if z != 1 else h_types.get(neighbors[k][0], h_default) for k, (_, z, _) in enumerate(lig.atoms)]
        ref_classes = [prm.atoms[t]["symbol"] for t in types]
        ref = np.array([prm.charges[t] for t in types])
        raw = np.array(lig.charges)
        agree = sum(1 for a, b in zip(lig.classes, ref_classes) if a == b)
        agree_total += agree
        mismatches = sorted({"%s->%s" % (b, a) for a, b in zip(lig.classes, ref_classes) if a != b})
        reference_q += list(ref); raw_q += list(raw)
        rows.append({"molecule": name, "atoms": len(lig.atoms), "classes_agree": agree, "mismatches": mismatches,
                     "mae_unscaled": float(np.mean(np.abs(raw - ref))), "mae_1.20": float(np.mean(np.abs(1.2 * raw - ref)))})
        print("%-18s %2d atoms  classes %2d/%2d %-12s MAE(1.00) %.3f  MAE(1.20) %.3f" % (
              name, len(lig.atoms), agree, len(lig.atoms), ",".join(mismatches) or "-", rows[-1]["mae_unscaled"], rows[-1]["mae_1.20"]))
    ref, raw = np.array(reference_q), np.array(raw_q)
    optimum = float(np.dot(raw, ref) / np.dot(raw, raw))
    print("\nclass agreement %d/%d   r = %.3f" % (agree_total, len(ref), np.corrcoef(raw, ref)[0, 1]))
    scales = {}
    for s in (1.0, optimum, 1.20, 1.27):
        q = s * raw
        scales["%.3f" % s] = {"mae": float(np.mean(np.abs(q - ref))), "rmsd": float(np.sqrt(np.mean((q - ref) ** 2)))}
        print("scale %.3f   MAE %.3f   RMSD %.3f" % (s, scales["%.3f" % s]["mae"], scales["%.3f" % s]["rmsd"]))
    if "--json" in sys.argv:
        json.dump({"molecules": rows, "scales": scales, "optimum": optimum},
                  open(sys.argv[sys.argv.index("--json") + 1], "w"), indent=1)


if __name__ == "__main__":
    main()
