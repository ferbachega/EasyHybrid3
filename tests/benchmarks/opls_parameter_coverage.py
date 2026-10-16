"""
Benchmark: coverage of the OPLS parameter database for small molecules.

Runs the full preparation (src/pdynamo/opls/prep.py, ligand parametrization with
CM5) on the .mol2 files of selected Structure Library categories and counts, per
term type, how many bonded parameters came from the database and how many were
estimated or set to zero.

Needs: pDynamo3, obabel, xtb.
Usage: python3 tests/benchmarks/opls_parameter_coverage.py [category ...] [--json out.json]
Default categories: those of the 2026-10-04 run (144 molecules: 140 prepared;
bonds 98.7 % from the database, angles 94.5 %, dihedrals 91.7 %, impropers 49 %).
"""
import collections, contextlib, glob, io, json, os, shutil, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path += [os.environ.get("PDYNAMO3_HOME", os.path.expanduser("~/programs/pDynamo3")),
             os.path.join(REPO, "src"), os.path.join(REPO, "src", "graphics_engine", "src")]

import pBabel                                    # noqa: E402,F401  (pDynamo import order)
from pBabel import ImportSystem                  # noqa: E402
from pdynamo.opls import prep, mol2_import       # noqa: E402

LIBRARY = os.path.join(REPO, "src", "gui", "windows", "builder", "structures")
DEFAULT = ["08_classic_fused_medchem_systems", "10_carbonyl_scaffolds", "27_vitamins",
           "28_neurotransmitters_hormones", "15_nucleosides", "05_saturated_n_heterocycles"]


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out_json = sys.argv[sys.argv.index("--json") + 1] if "--json" in sys.argv else None
    if out_json in args: args.remove(out_json)
    files = []
    for category in (args or DEFAULT):
        files += sorted(glob.glob(os.path.join(LIBRARY, category, "*.mol2")))
    results = []
    for path in files:
        name, start = os.path.relpath(path, LIBRARY), time.time()
        work = tempfile.mkdtemp(prefix="opls_cov_")
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                system, _ = mol2_import.import_system(path, log=None)
                _, report = prep.prepare_opls_system(system, work_dir=work)
            estimated = collections.Counter()
            for note in report.get("estimated_terms", []):
                kind = "dihedral_ring" if "aromatic ring" in note else ("dihedral_zero" if "set to zero" in note else note.split()[0])
                estimated[kind] += 1
            results.append({"file": name, "ok": report["ok"], "atoms": report["atoms"],
                            "sources": {t: dict(c) for t, c in report.get("parameters", {}).get("sources", {}).items()},
                            "estimated": dict(estimated), "time": time.time() - start})
        except Exception as error:
            results.append({"file": name, "ok": False, "error": "%s: %s" % (type(error).__name__, str(error)[:160]),
                            "time": time.time() - start})
        finally:
            shutil.rmtree(work, ignore_errors=True)
        print("%-60s %s" % (name, "ok" if results[-1]["ok"] else results[-1].get("error", "not ok")), flush=True)
    ok = [r for r in results if r["ok"]]
    print("\nmolecules %d, prepared %d" % (len(results), len(ok)))
    per_term = collections.defaultdict(collections.Counter)
    for r in ok:
        for term, counts in r["sources"].items():
            per_term[term].update(counts)
    for term, counts in per_term.items():
        total = sum(counts.values())
        print("%-11s %5d  " % (term, total) + "  ".join("%s %d (%.1f%%)" % (s, n, 100.0 * n / total) for s, n in counts.most_common()))
    if out_json:
        json.dump(results, open(out_json, "w"), indent=1)


if __name__ == "__main__":
    main()
