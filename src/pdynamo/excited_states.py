"""
excited_states
==============

Vertical excited states of a molecule (gas phase) with pDynamo:

  MNDO-CI   QCModelMNDOCI: configuration interaction (singles, doubles,
            singles + doubles or full) in an active space of a semiempirical
            (AM1, PM3, PM6, RM1, MNDO) reference. State energies and spins.
  MOM       Maximum overlap method (Gilbert, Besley & Gill 2008): an excited
            determinant made by swapping an occupied and a virtual orbital
            (pSimulation.MakeExcitedStateByOrbitalSwapping) and re-converged
            (SetUpMOMExcitedStateCalculation); the transition dipole moment
            and the oscillator strength from pSimulation.TransitionDipoleMoment.
            Semiempirical or HF/DFT (pdynamo/esp_charges.py levels).

Runs in a separate process from the GUI: python -m pdynamo.excited_states job.json result.json
"""

import json
import math
import os
import time
import traceback

KJ_PER_EV = 96.48533212
HARTREE_KJ = 2625.4996394799
EV_NM = 1239.84198

SEMIEMPIRICAL = ("am1", "pm3", "pm6", "rm1", "mndo")
CI_METHODS = ("Singles", "Doubles", "SinglesDoubles", "Full")


def _system(Z, xyz, charge, multiplicity, occupancy=None, restricted=True):
    import pBabel  # noqa: F401
    from pMolecule import System
    from pMolecule.QCModel import ElectronicState, OccupancyType
    from pScientific.Geometry3 import Coordinates3
    s = System.FromAtoms([int(z) for z in Z])
    c = Coordinates3.WithExtent(len(Z))
    for i, (x, y, z) in enumerate(xyz):
        c[i, 0], c[i, 1], c[i, 2] = float(x), float(y), float(z)
    s.coordinates3 = c
    opts = dict(charge=int(charge), multiplicity=int(multiplicity), isSpinRestricted=restricted)
    if occupancy == "fractional":
        opts["occupancyType"] = OccupancyType.FractionalFixed
        opts["numberFractionalHOOs"] = 1
        opts["numberFractionalLUOs"] = 1 if int(multiplicity) == 1 else 0
    s.electronicState = ElectronicState.WithOptions(**opts)
    return s


def _qc_model(method, ci=None):
    from pMolecule.QCModel import DIISSCFConverger, QCModelDFT, QCModelMNDO, QCModelMNDOCI, CIMethod
    converger = DIISSCFConverger.WithOptions(densityTolerance=1.0e-8, maximumIterations=500)
    if method.lower() in SEMIEMPIRICAL:
        if ci is not None:
            return QCModelMNDOCI.WithOptions(converger=converger, hamiltonian=method.lower(),
                                             ciMethod=getattr(CIMethod, ci["method"]),
                                             activeElectrons=int(ci["electrons"]), activeOrbitals=int(ci["orbitals"]),
                                             multiplicity=int(ci.get("multiplicity", 1)),
                                             minimalMultiplicity=int(ci.get("minimal", 1)),
                                             numberOfStates=int(ci["states"]))
        return QCModelMNDO.WithOptions(converger=converger, hamiltonian=method.lower())
    from pdynamo.esp_charges import LEVELS
    functional, basis = LEVELS[method]
    return QCModelDFT.WithOptions(converger=converger, functional=functional, orbitalBasis=basis)


def mndo_ci(Z, xyz, charge=0, multiplicity=1, hamiltonian="am1", method="SinglesDoubles",
            electrons=4, orbitals=4, states=10, minimal=1):
    """CI states (energies relative to the lowest state, kJ/mol and eV) and spins."""
    # . a closed-shell reference with fractional HOMO/LUMO occupancy, as in pDynamo's examples
    s = _system(Z, xyz, charge, multiplicity, occupancy="fractional")
    s.DefineQCModel(_qc_model(hamiltonian, {"method": method, "electrons": electrons, "orbitals": orbitals,
                                            "states": states, "multiplicity": multiplicity, "minimal": minimal}))
    s.Energy(log=None)
    node = s.scratch.ci
    energies = [node.ciEnergies[i] for i in range(len(node.ciEnergies))]
    spins = [node.ciSpins[i] for i in range(len(node.ciSpins))]
    e0 = min(energies)
    out = []
    for e, x in zip(energies, spins):
        S = 0.5 * (-1.0 + math.sqrt(max(1.0 + 4.0 * x, 0.0)))
        de = (e - e0) / KJ_PER_EV
        out.append({"energy_kj": e - e0, "energy_ev": de, "nm": EV_NM / de if de > 1e-6 else None,
                    "S": S, "multiplicity": int(round(2 * S + 1))})
    return {"states": out, "method": method, "hamiltonian": hamiltonian,
            "active": [int(electrons), int(orbitals)], "absolute_kj": e0}


def mom(Z, xyz, swaps, charge=0, multiplicity=1, method="am1"):
    """Ground state + one MOM excited state per swap (o, v) = (HOMO-o, LUMO+v)."""
    import pBabel  # noqa: F401  (import order: pBabel before pMolecule)
    from pCore import Clone
    from pMolecule.QCModel import ElectronicState, OccupancyType
    from pSimulation import MakeExcitedStateByOrbitalSwapping, SetUpMOMExcitedStateCalculation, TransitionDipoleMoment
    s = _system(Z, xyz, charge, multiplicity, restricted=False)
    s.DefineQCModel(_qc_model(method))
    e_gs = s.Energy(log=None)
    gs = (Clone(s.scratch.orbitalsP), Clone(s.scratch.orbitalsQ))
    homo = gs[0].occupancyHandler.numberOccupied - 1
    occP = [gs[0].energies[i] for i in range(len(gs[0].energies))]
    es_state = ElectronicState.WithOptions(charge=int(charge), multiplicity=int(multiplicity),
                                           isSpinRestricted=False, occupancyType=OccupancyType.MOM)
    SetUpMOMExcitedStateCalculation(s, es_state)
    results = []
    for (o, v) in swaps:
        gs[0].orbitals.CopyTo(s.scratch.orbitalsP.orbitals)
        gs[1].orbitals.CopyTo(s.scratch.orbitalsQ.orbitals)
        entry = {"swap": [int(o), int(v)], "label": "HOMO{} -> LUMO{}".format("-%d" % o if o else "", "+%d" % v if v else "")}
        try:
            MakeExcitedStateByOrbitalSwapping(s, swapsP=[[int(o), int(v)]])
            e_es = s.Energy(log=None)
            es = (Clone(s.scratch.orbitalsP), Clone(s.scratch.orbitalsQ))
            de = e_es - e_gs
            overlap, tdm, f = TransitionDipoleMoment(s, de / HARTREE_KJ, gs, es)
            ev = de / KJ_PER_EV
            entry.update({"energy_kj": de, "energy_ev": ev, "nm": EV_NM / ev if ev > 1e-6 else None,
                          "f": float(f), "tdm_debye": [float(tdm[k]) for k in range(3)],
                          "tdm_norm": math.sqrt(sum(float(tdm[k]) ** 2 for k in range(3))),
                          "overlap": float(overlap)})
            if abs(overlap) > 0.9 or ev < 0.05:
                entry["collapsed"] = True
        except Exception as error:
            entry["error"] = str(error)
        results.append(entry)
    return {"ground_kj": e_gs, "homo_index": homo, "orbital_energies": occP, "excited": results, "method": method}


def run_job_file(job_path, result_path):
    with open(job_path) as handle:
        job = json.load(handle)
    t0 = time.time()
    try:
        if job["kind"] == "ci":
            result = mndo_ci(job["Z"], job["xyz"], job["charge"], job["multiplicity"], job["hamiltonian"],
                             job["method"], job["electrons"], job["orbitals"], job["states"], job.get("minimal", 1))
        else:
            result = mom(job["Z"], job["xyz"], job["swaps"], job["charge"], job["multiplicity"], job["method"])
        result["seconds"] = time.time() - t0
    except Exception as error:
        result = {"error": "{}\n\n{}".format(error, traceback.format_exc(limit=4))}
    with open(result_path + ".tmp", "w") as handle:
        json.dump(result, handle)
    os.replace(result_path + ".tmp", result_path)


if __name__ == "__main__":
    import sys
    run_job_file(sys.argv[1], sys.argv[2])
