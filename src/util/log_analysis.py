#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Job log analysis (2026-10-28): the "Plots" tabs of the log window (double
#  click on a job of the Process Manager) and the job_summary command of the
#  automation layer. GTK-free: reads the pDynamo log of a job (and the ORCA /
#  xTB files next to it) and returns
#
#      {"kind": "optimization" | "md" | "scan1d" | "scan2d" | "neb" | "umbrella"
#               | "normal_modes" | "refinement" | "energy" | None,
#       "summary": [(label, value text), ...],
#       "panels":  [panel, ...],          # the pDynamo log
#       "qc":      {"program": "ORCA" | "xTB", "summary": [...], "panels": [...]} | None}
#
#  A panel is what util/analysis_plot.AnalysisPlot draws:
#      {"title", "kind": "lines", "series": [{"x", "y", "label", "style"}], "xlabel", "ylabel", "hline"}
#      {"title", "kind": "heatmap", "matrix", "xlabel", "ylabel", "row_labels", "col_labels",
#       "value_label", "overlays"}
#  Energies in kJ/mol unless said otherwise.

import glob
import math
import os
import re

_NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
HARTREE = 2625.499639


def _floats(text):
    try:
        return [float(v) for v in text.split()]
    except ValueError:
        return None


def _series(x, y, label, style="line", **extra):
    s = {"x": list(x), "y": list(y), "label": label, "style": style}
    s.update(extra)
    return s


def _lines_panel(title, series, xlabel, ylabel, hline=None):
    return {"title": title, "kind": "lines", "series": series, "xlabel": xlabel, "ylabel": ylabel, "hline": hline}


def _fmt(v, unit="", digits=2):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "-"
    return ("{:.%df}" % digits).format(v) + (" " + unit if unit else "")


# ============================================================================ pDynamo log pieces
def minimizer_runs(lines):
    """ Every minimizer table: [{"iter", "f", "rms", "status"}], in order. """
    runs, i, n = [], 0, len(lines)
    while i < n:
        if re.match(r"\s*Iteration\s+Function\s+RMS Gradient", lines[i]):
            it, f, rms = [], [], []
            i += 2
            while i < n and not lines[i].startswith("-----"):
                parts = lines[i].split()
                vals = [p for p in parts if re.fullmatch(_NUM, p)]
                if len(vals) >= 3 and parts[0].isdigit():
                    it.append(int(parts[0]))
                    f.append(float(vals[1]))
                    rms.append(float(vals[2]))
                i += 1
            status = None
            for j in range(i, min(i + 6, n)):
                text = lines[j].strip()
                if text and not text.startswith("-"):
                    status = text
                    break
            if it:
                runs.append({"iter": it, "f": f, "rms": rms, "status": status})
        i += 1
    return runs


def energy_terms(lines):
    """ The last "Summary of Energy Terms" block: {term: value}. """
    start = None
    for i, line in enumerate(lines):
        if "Summary of Energy Terms" in line:
            start = i
    if start is None:
        return {}
    terms = {}
    for line in lines[start + 1:start + 40]:
        if line.startswith("-----") and terms:
            break
        for name, value in re.findall(r"([A-Za-z][A-Za-z0-9 /.\-()]*?)\s*=\s*(" + _NUM + r"|None)", line):
            if value != "None":
                terms[name.strip()] = float(value)
    return terms


def md_table(lines):
    """ The dynamics tables (Time, Total/Kinetic/Potential Energy, Temperature, [Pressure, Volume]) joined. """
    data, header = {}, None
    for i, line in enumerate(lines):
        if re.match(r"\s*Time\s+Total Energy", line):
            header = [h.strip() for h in re.split(r"\s{2,}", line.strip())]
            continue
        if header is None:
            continue
        if line.startswith("-----"):
            continue
        vals = _floats(line)
        if vals and len(vals) == len(header):
            for h, v in zip(header, vals):
                data.setdefault(h, []).append(v)
        elif line.strip() and not vals:
            header = None
    return data


# ============================================================================ job kinds
def _analyze_optimization(lines, out):
    runs = minimizer_runs(lines)
    if not runs:
        return False
    r = runs[-1]
    out["kind"] = "optimization"
    out["summary"] += [("status", r["status"] or "-"), ("iterations", str(r["iter"][-1])),
                       ("final energy", _fmt(r["f"][-1], "kJ/mol", 4)),
                       ("change from start", _fmt(r["f"][-1] - r["f"][0], "kJ/mol")),
                       ("final RMS gradient", _fmt(r["rms"][-1], "kJ/mol/A", 4))]
    out["panels"].append(_lines_panel("Energy along the optimization",
                                      [_series(r["iter"], r["f"], "energy", "line")], "iteration", "energy (kJ/mol)"))
    out["panels"].append(_lines_panel("RMS gradient (log10)",
                                      [_series(r["iter"], [math.log10(max(v, 1e-12)) for v in r["rms"]], "log10 RMS gradient")],
                                      "iteration", "log10 (kJ/mol/A)"))
    return True


def _analyze_md(lines, out):
    d = md_table(lines)
    if "Time" not in d or len(d["Time"]) < 2:
        return False
    out["kind"] = "md"
    t = d["Time"]
    stats = []
    for key in ("Temperature", "Total Energy", "Potential Energy", "Pressure", "Volume"):
        if key in d:
            v = d[key][len(d[key]) // 5:] or d[key]           # . the first fifth taken as relaxation
            mean = sum(v) / len(v)
            std = (sum((x - mean) ** 2 for x in v) / len(v)) ** 0.5
            stats.append((key.lower(), "{:.2f} ± {:.2f}".format(mean, std)))
    drift = None
    if "Total Energy" in d and t[-1] > t[0]:
        drift = (d["Total Energy"][-1] - d["Total Energy"][0]) / (t[-1] - t[0])
    out["summary"] += [("time", "{:.3f} ps, {} points".format(t[-1], len(t)))] + stats + \
                      [("total energy drift", _fmt(drift, "kJ/mol/ps"))] + \
                      [("note", "means and std over the last 4/5 of the run")]
    out["panels"].append(_lines_panel("Energies", [_series(t, d[k], k) for k in ("Total Energy", "Potential Energy", "Kinetic Energy") if k in d],
                                      "time (ps)", "kJ/mol"))
    if "Temperature" in d:
        out["panels"].append(_lines_panel("Temperature", [_series(t, d["Temperature"], "T")], "time (ps)", "K"))
    for key, unit in (("Pressure", "atm"), ("Volume", "A^3")):
        if key in d:
            out["panels"].append(_lines_panel(key, [_series(t, d[key], key)], "time (ps)", unit))
    return True


def _profile_panels(out, rc, z, xlabel, title):
    from automation.api import _profile_1d
    prof = _profile_1d(rc, z)
    out["summary"] += [("points", str(len(z))),
                       ("highest point", "{} at {}".format(_fmt(prof["highest"]["energy"], "kJ/mol"), prof["highest"]["rc"]))]
    for m in prof.get("local_maxima", []):
        out["summary"].append(("barrier (local max)", "{} at {}".format(_fmt(m["energy"], "kJ/mol"), m["rc"])))
    for m in prof.get("local_minima", []):
        out["summary"].append(("local minimum", "{} at {}".format(_fmt(m["energy"], "kJ/mol"), m["rc"])))
    out["summary"].append(("first → last", _fmt(prof.get("first_to_last"), "kJ/mol")))
    out["panels"].append(_lines_panel(title, [_series(rc, z, "energy", "line"), _series(rc, z, "points", "points")],
                                      xlabel, "relative energy (kJ/mol)"))


def _analyze_scan(logfile, lines, out):
    from pdynamo.LogFileWriter import LogFileReader
    try:
        data = LogFileReader(logfile).get_data()
    except Exception:
        data = None
    if not data:
        return False
    import numpy as np
    Z = np.asarray(data["Z"], dtype=float)
    Z = Z - np.nanmin(Z)
    if data.get("type") == "plot2D":
        from automation.api import _surface_2d
        Z = Z.T
        rc1, rc2 = np.asarray(data["RC1"], dtype=float).T, np.asarray(data["RC2"], dtype=float).T
        out["kind"] = "scan2d" if "SCAN2D" in "".join(lines[:4000]) else "refinement"
        surf = _surface_2d(Z, lambda i, j: rc1[i, j], lambda i, j: rc2[i, j])
        out["summary"] += [("grid", "{} x {} (RC1 x RC2)".format(*Z.shape)),
                           ("lowest", "cell {} ({}, {})".format(surf["lowest"]["cell"], surf["lowest"]["rc1"], surf["lowest"]["rc2"])),
                           ("highest", "{} kJ/mol".format(surf["highest"]["energy"]))]
        overlays = []
        if surf.get("path"):
            p = surf["path"]
            out["summary"] += [("path barrier (minimum to minimum)", _fmt(p["barrier"], "kJ/mol")),
                               ("path top", "cell {} ({}, {})".format(p["top"]["cell"], p["top"]["rc1"], p["top"]["rc2"])),
                               ("reaction energy", _fmt(p["reaction_energy"], "kJ/mol"))]
            overlays = [{"cells": [(i, j) for i, j in p["cells"]], "style": "path", "color": (1, 1, 1), "width": 2},
                        {"cells": [tuple(p["from"]["cell"]), tuple(p["to"]["cell"])], "style": "dots", "color": (0, 0, 0), "radius": 4},
                        {"cells": [tuple(p["top"]["cell"])], "style": "dots", "color": (0.85, 0.1, 0.1), "radius": 5}]
        out["panels"].append({"title": "Energy surface: lowest-barrier path (white), its ends (black), its top (red)", "kind": "heatmap",
                              "matrix": [[float(v) for v in row] for row in Z], "xlabel": "RC2", "ylabel": "RC1",
                              "row_labels": ["{:.2f}".format(v) for v in rc1[:, 0]],
                              "col_labels": ["{:.2f}".format(v) for v in rc2[0, :]],
                              "value_label": "kJ/mol", "overlays": overlays})
        return True
    rc = [float(v) for v in data["RC1"]]
    out["kind"] = "scan1d" if "EasyHybrid-SCAN" in "".join(lines[:4000]) else (
        "neb" if any("Path Summary" in l for l in lines) else "refinement")
    if out["kind"] == "neb":
        return False
    _profile_panels(out, rc, list(Z), "reaction coordinate", "Energy profile")
    if out["kind"] == "scan1d":
        runs = minimizer_runs(lines)
        if len(runs) >= len(rc) - 1 and runs:
            runs = runs[-len(rc):] if len(runs) > len(rc) else runs
            xs = rc[-len(runs):]
            out["panels"].append(_lines_panel("Optimization of each point",
                                              [_series(xs, [r["iter"][-1] for r in runs], "iterations", "bars")],
                                              "reaction coordinate", "iterations"))
            bad = [x for x, r in zip(xs, runs) if r["status"] and "converged" not in r["status"].lower()]
            out["summary"].append(("points not converged", "{} of {}{}".format(
                len(bad), len(runs), (": RC " + ", ".join("{:.2f}".format(b) for b in bad[:6])) if bad else "")))
    return True


def _analyze_neb(lines, out):
    start = None
    for i, line in enumerate(lines):
        if "Path Summary" in line:
            start = i
    if start is None:
        return False
    rows = []
    for line in lines[start + 4:]:
        vals = _floats(line)
        if not vals or len(vals) < 9:
            if rows:
                break
            continue
        rows.append(vals)
    if not rows:
        return False
    out["kind"] = "neb"
    image = [int(r[0]) for r in rows]
    rel = [r[2] for r in rows]
    rms = [r[3] for r in rows]
    arc = [r[-2] for r in rows]
    top = max(range(len(rel)), key=lambda k: rel[k])
    out["summary"] += [("images", str(len(rows))), ("highest image", "{} ({})".format(image[top], _fmt(rel[top], "kJ/mol"))),
                       ("barrier from the first image", _fmt(rel[top] - rel[0], "kJ/mol")),
                       ("last − first", _fmt(rel[-1] - rel[0], "kJ/mol")),
                       ("largest RMS gradient (inner images)", _fmt(max(rms[1:-1] or [0.0]), "kJ/mol/A", 3))]
    out["panels"].append(_lines_panel("Energy along the path", [_series(arc, rel, "energy"), _series(arc, rel, "images", "points")],
                                      "arc length", "relative energy (kJ/mol)"))
    out["panels"].append(_lines_panel("RMS gradient of each image", [_series(image, rms, "RMS gradient", "bars")],
                                      "image", "kJ/mol/A"))
    return True


def _analyze_normal_modes(lines, out):
    def block(title):
        for i, line in enumerate(lines):
            if title in line:
                vals = []
                for l2 in lines[i + 2:]:
                    if l2.startswith("-----"):
                        break
                    v = _floats(l2)
                    if v is None:
                        break
                    vals += v
                return vals
        return None
    freqs = block("Harmonic Frequencies")
    if not freqs:
        return False
    out["kind"] = "normal_modes"
    intens = block("Infrared Intensities")
    imag = [f for f in freqs if f < -20.0]
    real = [f for f in freqs if f > 20.0]
    out["summary"] += [("frequencies", str(len(freqs))),
                       ("imaginary (< −20 cm⁻¹)", "{}{}".format(len(imag), (": " + ", ".join("{:.1f}".format(f) for f in imag[:5])) if imag else "")),
                       ("lowest real", _fmt(min(real) if real else None, "cm⁻¹", 1))]
    y = intens if intens and len(intens) == len(freqs) else [1.0] * len(freqs)
    out["panels"].append(_lines_panel("Frequencies" + (" (IR intensity)" if y is intens else ""),
                                      [_series(freqs, y, "modes", "bars")], "wavenumber (cm⁻¹)",
                                      "IR intensity" if y is intens else "mode"))
    return True


def _analyze_umbrella(logfile, out):
    folder = os.path.dirname(logfile)
    if not os.path.isdir(os.path.join(folder, "data_collection")):
        return False
    from pdynamo.us_data import USRun
    run = USRun(folder)
    out["kind"] = "umbrella"
    stats = run.stats()
    keys = [k for k in run.keys if (stats.get(k) or {}).get("n")]
    dev = [stats[k]["deviation"] for k in keys]
    out["summary"] += [("windows", "{} ({}D)".format(len(run.keys), run.dims)), ("temperature", _fmt(run.temperature, "K", 1))]
    if run.dims == 1:
        pairs = run.overlap_1d()[0]
        ov = [o for _a, _b, o in pairs]
        low = [o for o in ov if o < 0.03]
        out["summary"] += [("lowest overlap of neighbours", _fmt(min(ov) if ov else None, "", 3)),
                           ("pairs below 0.03", "{} of {}".format(len(low), len(ov))),
                           ("windows > 3 σ from target", str(sum(1 for d in dev if d > 3.0))),
                           ("verdict", "fail" if low or any(d > 3.0 for d in dev) else "OK"),
                           ("more", "Analysis > Umbrella Sampling Browser")]
        out["panels"].append(_lines_panel("Overlap of neighbouring windows (needs ≥ 0.03)",
                                          [_series(list(range(len(ov))), ov, "overlap", "bars")], "pair", "overlap", hline=0.03))
        out["panels"].append(_lines_panel("|mean − target| / restraint width (≤ 3)",
                                          [_series([run.target(k)[0] for k in keys], dev, "deviation", "bars")],
                                          "window target", "σ", hline=3.0))
    return True


def _analyze_energy(lines, out):
    terms = energy_terms(lines)
    if not terms:
        return False
    out["kind"] = out["kind"] or "energy"
    total = terms.pop("Potential Energy", None)
    out["summary"].append(("potential energy", _fmt(total, "kJ/mol", 4)))
    names = list(terms)
    out["panels"].append(_lines_panel("Energy terms (last evaluation)",
                                      [_series(list(range(len(names))), [terms[k] for k in names], "terms", "bars")],
                                      "", "kJ/mol"))
    out["panels"][-1]["xlabels"] = names
    return True


# ============================================================================ QC programs
def _orca_files(folder):
    pats = ("*.orca.log", "*.orca.out", "orcaJob.log", "*orca*.log")
    return sorted({p for pat in pats for p in glob.glob(os.path.join(folder, pat)) if not os.path.basename(p).startswith("error_")})


def _xtb_files(folder):
    pats = ("*.xtb.out", "*.xtb.log", "XTBJob*.out", "XTBJob*.log")
    return sorted({p for pat in pats for p in glob.glob(os.path.join(folder, pat))})


def _natural_key(path):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", os.path.basename(path))]


def analyze_orca(paths):
    energies, last_scf, labels, gaps = [], None, [], []
    for path in sorted(paths, key=_natural_key):
        try:
            text = open(path, errors="replace").read()
        except OSError:
            continue
        for m in re.finditer(r"FINAL SINGLE POINT ENERGY\s+(" + _NUM + ")", text):
            energies.append(float(m.group(1)))
            labels.append(os.path.basename(path))
        scf = []
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if re.match(r"\s*(Iteration|ITER)\s+Energy", line):
                scf = []
                for l2 in lines[i + 1:]:
                    m = re.match(r"\s*(\d+)\s+(" + _NUM + r")\s+(" + _NUM + ")", l2)
                    if m:
                        scf.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
                    elif "SCF CONVERGED" in l2 or "SUCCESS" in l2 or "Total Energy" in l2:
                        break
        if scf:
            last_scf = scf
        occ = re.findall(r"^\s*\d+\s+(\d\.\d+)\s+(" + _NUM + r")\s+(" + _NUM + ")", text, re.M)
        if occ:
            homo = max((float(e) for o, _h, e in occ if float(o) > 0.5), default=None)
            lumo = min((float(e) for o, _h, e in occ if float(o) < 0.5), default=None)
            if homo is not None and lumo is not None:
                gaps.append(lumo - homo)
    if not energies and not last_scf:
        return None
    out = {"program": "ORCA", "summary": [("ORCA calculations found", str(len(energies)))], "panels": []}
    if energies:
        rel = [(e - min(energies)) * HARTREE for e in energies]
        out["summary"] += [("last energy", "{:.8f} Eh".format(energies[-1])),
                           ("span of the energies", _fmt(max(rel), "kJ/mol"))]
        out["panels"].append(_lines_panel("ORCA energy of each calculation (relative)",
                                          [_series(list(range(len(rel))), rel, "energy", "line"),
                                           _series(list(range(len(rel))), rel, "calculations", "points")],
                                          "calculation (file order)", "kJ/mol"))
    if last_scf:
        it = [s[0] for s in last_scf]
        de = [math.log10(max(abs(s[2]), 1e-14)) for s in last_scf[1:]]
        out["summary"].append(("last SCF", "{} cycles".format(len(last_scf))))
        if de:
            out["panels"].append(_lines_panel("Last SCF: log10 |ΔE| per cycle", [_series(it[1:], de, "log10 |ΔE|")],
                                              "cycle", "log10 (Eh)"))
    if gaps:
        out["summary"].append(("last HOMO–LUMO gap", _fmt(gaps[-1], "eV")))
    return out


def analyze_xtb(paths):
    energies, gaps = [], []
    for path in sorted(paths, key=_natural_key):
        try:
            text = open(path, errors="replace").read()
        except OSError:
            continue
        e = re.findall(r"TOTAL ENERGY\s+(" + _NUM + r")\s+Eh", text)
        g = re.findall(r"HOMO-LUMO gap\s+(" + _NUM + r")\s+eV", text, re.I)
        if e:
            energies.append(float(e[-1]))
            gaps.append(float(g[-1]) if g else float("nan"))
    if not energies:
        return None
    rel = [(e - min(energies)) * HARTREE for e in energies]
    xs = list(range(len(energies)))
    out = {"program": "xTB", "summary": [("xTB calculations found", str(len(energies))),
                                          ("last energy", "{:.8f} Eh".format(energies[-1])),
                                          ("span of the energies", _fmt(max(rel), "kJ/mol")),
                                          ("HOMO–LUMO gap", "{:.2f} – {:.2f} eV".format(
                                              min(g for g in gaps if g == g), max(g for g in gaps if g == g))
                                           if any(g == g for g in gaps) else "-")],
           "panels": [_lines_panel("xTB energy of each calculation (relative)",
                                   [_series(xs, rel, "energy"), _series(xs, rel, "calculations", "points")],
                                   "calculation (file order)", "kJ/mol")]}
    if any(g == g for g in gaps):
        out["panels"].append(_lines_panel("HOMO–LUMO gap", [_series(xs, gaps, "gap", "points")],
                                          "calculation (file order)", "eV"))
    return out


# ============================================================================ entry point
def analyze(logfile, simulation_type=None):
    """ Analysis of a job log (see the module docstring). Never raises for an
        unknown or partial log: kind None and empty lists then. """
    out = {"kind": None, "summary": [], "panels": [], "qc": None, "logfile": logfile}
    try:
        lines = open(logfile, errors="replace").read().splitlines()
    except OSError as error:
        out["summary"].append(("error", str(error)))
        return out
    st = simulation_type or ""
    folder = os.path.dirname(logfile)
    try:
        if st == "Umbrella_Sampling" or (not st and os.path.isdir(os.path.join(folder, "data_collection"))):
            _analyze_umbrella(logfile, out)
        elif st == "Molecular_Dynamics" or (not st and md_table(lines)):
            _analyze_md(lines, out)
        elif st == "Nudged_Elastic_Band" or any("Path Summary" in l for l in lines[-400:]):
            _analyze_neb(lines, out)
        elif st == "Normal_Modes" or (not st and any("Harmonic Frequencies" in l for l in lines)):
            _analyze_normal_modes(lines, out)
        elif st in ("Relaxed_Surface_Scan", "Advanced_Relaxed_Surface_Scan", "Energy_Refinement") or \
                any(l.startswith("TYPE") for l in lines[:2000]):
            _analyze_scan(logfile, lines, out)
        if out["kind"] is None:
            if not _analyze_optimization(lines, out):
                _analyze_energy(lines, out)
        elif out["kind"] in ("optimization", "neb", "normal_modes"):
            pass
    except Exception as error:                       # . a partial log of a running job, say
        out["summary"].append(("could not analyze", "{}: {}".format(type(error).__name__, error)))
    if out["kind"] in ("energy", "optimization") and not any(p["title"].startswith("Energy terms") for p in out["panels"]):
        _analyze_energy(lines, out)
    # . QC programs next to the log (one file per scan point / calculation)
    try:
        orca = _orca_files(folder)
        xtb = _xtb_files(folder)
        out["qc"] = analyze_orca(orca) if orca else (analyze_xtb(xtb) if xtb else None)
    except Exception as error:
        out["qc"] = {"program": "QC", "summary": [("could not analyze", str(error))], "panels": []}
    return out
