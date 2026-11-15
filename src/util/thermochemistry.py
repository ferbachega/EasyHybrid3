"""
thermochemistry
===============

Ideal-gas, rigid-rotor, harmonic-oscillator (RRHO) thermochemistry from
vibrational frequencies, as in J. W. Ochterski, "Thermochemistry in
Gaussian" (2000) and McQuarrie, Statistical Mechanics.

pSimulation.ThermodynamicsRRHO_SystemGeometry is not used: its final
"ideal gas" correction subtracts R (ln N_A - 1) from the entropy although
its translational partition function is already per molecule (v = kT/p),
which gives S(water, 298 K) = -267 J/(mol K) instead of +189, and its
enthalpy is never set.

Optionally the quasi-RRHO treatment of Grimme (Chem. Eur. J. 18, 9955,
2012) replaces the vibrational entropy of low modes by a free-rotor entropy
(damping function with w0 = 100 cm^-1, alpha = 4).

    thermo = rrho(frequencies_cm, masses_amu, coords_A, temperature=298.15,
                  pressure_atm=1.0, symmetry_number=1, multiplicity=1)
"""

import math

import numpy as np

H = 6.62607015e-34          # . J s
KB = 1.380649e-23           # . J/K
NA = 6.02214076e23
C = 2.99792458e10           # . cm/s
AMU = 1.66053906660e-27     # . kg
ATM = 101325.0              # . Pa
R = KB * NA                 # . J/(mol K)

# . rotational symmetry numbers of the point groups
SYMMETRY_NUMBERS = {"C1": 1, "Ci": 1, "Cs": 1, "C*v": 1, "Cinfv": 1, "D*h": 2, "Dinfh": 2,
                    "C2": 2, "C2v": 2, "C2h": 2, "C3": 3, "C3v": 3, "C3h": 3, "C4": 4, "C4v": 4, "C4h": 4,
                    "C5": 5, "C5v": 5, "C5h": 5, "C6": 6, "C6v": 6, "C6h": 6,
                    "D2": 4, "D2d": 4, "D2h": 4, "D3": 6, "D3d": 6, "D3h": 6, "D4": 8, "D4d": 8, "D4h": 8,
                    "D5": 10, "D5d": 10, "D5h": 10, "D6": 12, "D6d": 12, "D6h": 12,
                    "S4": 2, "S6": 3, "S8": 4, "T": 12, "Td": 12, "Th": 12, "O": 24, "Oh": 24, "I": 60, "Ih": 60}


def symmetry_number(point_group):
    return SYMMETRY_NUMBERS.get(str(point_group).strip(), 1)


def principal_moments(masses, coords):
    """Principal moments of inertia (amu A^2), ascending."""
    m = np.asarray(masses, dtype=float)
    x = np.asarray(coords, dtype=float)
    x = x - (m[:, None] * x).sum(axis=0) / m.sum()
    I = np.zeros((3, 3))
    for mi, (a, b, c) in zip(m, x):
        I += mi * np.array([[b * b + c * c, -a * b, -a * c], [-a * b, a * a + c * c, -b * c], [-a * c, -b * c, a * a + b * b]])
    return np.sort(np.linalg.eigvalsh(I))


def rrho(frequencies_cm, masses, coords, temperature=298.15, pressure_atm=1.0, symmetry_number=1,
         multiplicity=1, n_rt=None, quasi_rrho=False, w0=100.0, alpha=4.0, linear_tolerance=1e-3):
    """Thermochemistry (kJ/mol and J/(mol K)) of one molecule in the ideal gas.

    frequencies_cm: all 3N frequencies or only the vibrational ones; the
    n_rt lowest |frequencies| (6, or 5 for a linear molecule) are taken as
    rotations/translations when 3N values are given. Imaginary frequencies
    (negative values) among the vibrations are reported and left out.
    """
    T = float(temperature)
    m = np.asarray(masses, dtype=float)
    n_atoms = len(m)
    moments = principal_moments(m, coords) if n_atoms > 1 else np.zeros(3)
    linear = n_atoms == 2 or (n_atoms > 2 and moments[0] < linear_tolerance * max(moments[2], 1e-12))
    freqs = np.asarray(frequencies_cm, dtype=float)
    if n_rt is None:
        n_rt = 0 if n_atoms == 1 else (5 if linear else 6)
    if len(freqs) == 3 * n_atoms:
        order = np.argsort(np.abs(freqs))
        vib = np.sort(freqs[order[n_rt:]])
    else:
        vib = np.sort(freqs)
    imaginary = [float(f) for f in vib if f < 0]
    vib = vib[vib > 0]

    RT = R * T / 1000.0                 # . kJ/mol
    out = {"temperature": T, "pressure_atm": pressure_atm, "symmetry_number": symmetry_number,
           "linear": bool(linear), "imaginary": imaginary, "n_vib": int(len(vib))}

    # . translation
    mass = m.sum() * AMU
    v = KB * T / (pressure_atm * ATM)
    q_t = (2.0 * math.pi * mass * KB * T / H ** 2) ** 1.5 * v
    trans = {"S": R * (math.log(q_t) + 2.5), "E": 1.5 * RT, "Cv": 1.5 * R}
    # . electronic (one state of the given multiplicity)
    elec = {"S": R * math.log(multiplicity), "E": 0.0, "Cv": 0.0}
    # . rotation
    if n_atoms == 1:
        rot = {"S": 0.0, "E": 0.0, "Cv": 0.0}
    else:
        I = moments * AMU * 1e-20       # . kg m^2
        if linear:
            theta = H ** 2 / (8 * math.pi ** 2 * I[2] * KB)
            q_r = T / (symmetry_number * theta)
            rot = {"S": R * (math.log(q_r) + 1.0), "E": RT, "Cv": R}
        else:
            thetas = H ** 2 / (8 * math.pi ** 2 * I * KB)
            q_r = math.sqrt(math.pi) / symmetry_number * T ** 1.5 / math.sqrt(thetas.prod())
            rot = {"S": R * (math.log(q_r) + 1.5), "E": 1.5 * RT, "Cv": 1.5 * R}
    # . vibration (harmonic; energies from the bottom of the well)
    theta_v = H * C * vib / KB
    x = theta_v / T
    zpe = 0.5 * R * theta_v.sum() / 1000.0
    with np.errstate(over="ignore"):
        ex = np.exp(-x)
        e_vib = R * (theta_v * (0.5 + ex / (1 - ex))).sum() / 1000.0
        s_modes = R * (x * ex / (1 - ex) - np.log(1 - ex))
        cv = R * (x ** 2 * ex / (1 - ex) ** 2).sum()
    s_vib = s_modes.sum()
    if quasi_rrho and len(vib):
        # . Grimme: free-rotor entropy for low modes, interpolated by a damping function
        Bav = 1e-44                         # . kg m^2, average molecular moment
        mu = H / (8 * math.pi ** 2 * C * vib)   # . kg m^2 (cm in C cancels with cm^-1)
        mu_eff = mu * Bav / (mu + Bav)
        s_rot = R * (0.5 + np.log(np.sqrt(8 * math.pi ** 3 * mu_eff * KB * T / H ** 2)))
        wdamp = 1.0 / (1.0 + (w0 / vib) ** alpha)
        s_vib = (wdamp * s_modes + (1 - wdamp) * s_rot).sum()
        out["quasi_rrho"] = True
    vibd = {"S": s_vib, "E": e_vib, "Cv": cv}
    out["zpe"] = zpe

    parts = {"electronic": elec, "translational": trans, "rotational": rot, "vibrational": vibd}
    total = {k: sum(p[k] for p in parts.values()) for k in ("S", "E", "Cv")}
    total["H"] = total["E"] + RT
    total["G"] = total["H"] - T * total["S"] / 1000.0
    total["Cp"] = total["Cv"] + R
    out["parts"] = parts
    out["total"] = total
    out["moments"] = [float(v) for v in moments]
    out["per_mode_S"] = [float(s) for s in s_modes]
    out["vibrations"] = [float(f) for f in vib]
    return out


def format_table(t):
    """Text table of a rrho() result."""
    T = t["temperature"]
    lines = ["T = {:.2f} K, p = {:.3f} atm, symmetry number {}, {} vibrations{}{}".format(
        T, t["pressure_atm"], t["symmetry_number"], t["n_vib"], ", linear" if t["linear"] else "",
        ", quasi-RRHO (Grimme)" if t.get("quasi_rrho") else "")]
    if t["imaginary"]:
        lines.append("WARNING: {} imaginary frequencies left out: {}".format(
            len(t["imaginary"]), ", ".join("{:.1f}i".format(-f) for f in t["imaginary"])))
    lines.append("")
    lines.append("{:<15s} {:>14s} {:>14s} {:>14s}".format("", "E (kJ/mol)", "S (J/mol/K)", "Cv (J/mol/K)"))
    for name, p in t["parts"].items():
        lines.append("{:<15s} {:>14.3f} {:>14.3f} {:>14.3f}".format(name, p["E"], p["S"], p["Cv"]))
    tot = t["total"]
    lines.append("{:<15s} {:>14.3f} {:>14.3f} {:>14.3f}".format("total", tot["E"], tot["S"], tot["Cv"]))
    lines.append("")
    lines.append("Zero-point energy        ZPE = {:12.3f} kJ/mol".format(t["zpe"]))
    lines.append("Thermal energy           E   = {:12.3f} kJ/mol  (includes ZPE)".format(tot["E"]))
    lines.append("Thermal enthalpy         H   = {:12.3f} kJ/mol  (E + RT)".format(tot["H"]))
    lines.append("Entropy                  S   = {:12.3f} J/(mol K)".format(tot["S"]))
    lines.append("Thermal free energy      G   = {:12.3f} kJ/mol  (H - TS)".format(tot["G"]))
    lines.append("Heat capacity            Cp  = {:12.3f} J/(mol K)".format(tot["Cp"]))
    lines.append("")
    lines.append("Add the electronic energy to E, H and G for absolute values; G and H are the")
    lines.append("thermal corrections relative to the bottom of the well.")
    return "\n".join(lines)
